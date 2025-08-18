#ifndef MESH_H
#define MESH_H

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Macros.H"

#ifdef USE_OMP
    #include <omp.h>
#endif

enum BoundaryCondition
{
    TRANSMISSIVE = 0,
    REFLECTIVE = 1,
    PERIODIC = 2
};

class Geometry
{
    public:
        Geometry(const std::array<REAL, GRIDDIM>& lo, 
                 const std::array<REAL, GRIDDIM>& hi, 
                 const std::array<int, GRIDDIM>& res);
        const std::array<REAL, GRIDDIM>& getLo() const;
        const std::array<REAL, GRIDDIM>& getHi() const;
        const std::array<int, GRIDDIM>& getRes() const;
        const std::array<REAL, GRIDDIM>& getDx() const;
        void getPos(std::array<REAL, GRIDDIM>& pos, const std::array<int, GRIDDIM>& idx) const;
        void getIdx(std::array<REAL, GRIDDIM>& idx, const std::array<REAL, GRIDDIM>& pos) const;
        
    private:
        const std::array<REAL, GRIDDIM> m_lo;
        const std::array<REAL, GRIDDIM> m_hi;
        const std::array<int, GRIDDIM> m_res;
        const std::array<REAL, GRIDDIM> m_dx;
};

template<int NVARS>
class DataArray
{
    public:
        DataArray(const std::array<int, GRIDDIM>& res, const int nGhost=0) : m_res(res), m_nGhost(nGhost)
        {
            #if GRIDDIM == 1
                m_data.resize(res[0] + 2 * nGhost);
            #elif GRIDDIM == 2
                m_data.resize(res[0] + 2 * nGhost);
                for(int i = 0; i < res[0] + 2 * nGhost; ++i)
                {
                    m_data[i].resize(res[1] + 2 * nGhost);
                }
            #else  // GRIDDIM == 3
                m_data.resize(res[0] + 2 * nGhost);
                for(int i = 0; i < res[0] + 2 * nGhost; ++i)
                {
                    m_data[i].resize(res[1] + 2 * nGhost);
                    for(int j = 0; j < res[1] + 2 * nGhost; ++j)
                    {
                        m_data[i][j].resize(res[2] + 2 * nGhost);
                    }
                }
            #endif
        }

        const std::array<int, GRIDDIM>& getRes() const 
        {
            return m_res;
        }

        int getNGhost() const 
        {
            return m_nGhost;
        }

        std::array<REAL, NVARS>& operator()(GRIDDIM_DECL(const int i, const int j, const int k))
        {
            #if GRIDDIM == 1
                return m_data[i + m_nGhost];
            #elif GRIDDIM == 2
                return m_data[i + m_nGhost][j + m_nGhost];
            #else  // GRIDDIM == 3
                return m_data[i + m_nGhost][j + m_nGhost][k + m_nGhost];
            #endif
        }

        const std::array<REAL, NVARS>& operator()(GRIDDIM_DECL(const int i, const int j, const int k)) const
        {
            #if GRIDDIM == 1
                return m_data[i + m_nGhost];
            #elif GRIDDIM == 2
                return m_data[i + m_nGhost][j + m_nGhost];
            #else  // GRIDDIM == 3
                return m_data[i + m_nGhost][j + m_nGhost][k + m_nGhost];
            #endif
        }

        std::array<REAL, NVARS> operator()(GRIDDIM_DECL(const REAL i, const REAL j, const REAL k)) const
        {
            return (*this)({GRIDDIM_DECL(i, j, k)});
        }

        std::array<REAL, NVARS>& operator()(const std::array<int, GRIDDIM>& idx)
        {
            return (*this)(GRIDDIM_DECL(idx[0], idx[1], idx[2]));
        }

        const std::array<REAL, NVARS>& operator()(const std::array<int, GRIDDIM>& idx) const
        {
            return (*this)(GRIDDIM_DECL(idx[0], idx[1], idx[2]));
        }

        std::array<REAL, NVARS> operator()(const std::array<REAL, GRIDDIM>& idx) const
        {
            std::array<int, GRIDDIM> idxLo;
            std::array<std::array<REAL, 2>, GRIDDIM> fac;
            for(int d = 0; d < GRIDDIM; ++d)
            {
                idxLo[d] = static_cast<int>(std::floor(idx[d]));
                fac[d][1] = idx[d] - idxLo[d];
                fac[d][0] = 1.0 - fac[d][1];
            }

            std::array<REAL, NVARS> interpData;
            for(int v = 0; v < NVARS; ++v)
            {
                interpData[v] = 0.0;
            }
            
            for(int iOffset = 0; iOffset < 2; ++iOffset)
            {
                #if GRIDDIM >= 2
                for(int jOffset = 0; jOffset < 2; ++jOffset)
                #endif
                {
                    #if GRIDDIM == 3
                    for(int kOffset = 0; kOffset < 2; ++kOffset)
                    #endif
                    {
                        const std::array<REAL, NVARS> data = (*this)(GRIDDIM_DECL(idxLo[0] + iOffset, 
                                                                                  idxLo[1] + jOffset, 
                                                                                  idxLo[2] + kOffset));
                        for(int v = 0; v < NVARS; ++v)
                        {
                            interpData[v] += data[v] * GRIDDIM_TERM(fac[0][iOffset], 
                                                                    * fac[1][jOffset], 
                                                                    * fac[2][kOffset]);
                        }
                    }
                }
            }

            return interpData;
        }

        void fillGhost(const std::array<std::array<BoundaryCondition, GRIDDIM>, 2>& bc, 
                       const std::vector<std::array<int, GRIDDIM>>& vecIdx)
        {
            std::array<int, GRIDDIM> idxSource, idxTarget;
            std::array<int, GRIDDIM> nIter;
            for(int d = 0; d < GRIDDIM; ++d)
            {
                nIter = m_res;
                nIter[d] = m_nGhost;
                for(int s = 0; s < 2; ++s)
                {
                    for(int i = 0; i < nIter[0]; ++i)
                    {
                        #if GRIDDIM >= 2
                        for(int j = 0; j < nIter[1]; ++j)
                        #endif
                        {
                            #if GRIDDIM == 3
                            for(int k = 0; k < nIter[2]; ++k)
                            #endif
                            {
                                idxTarget = {GRIDDIM_DECL(i, j, k)};
                                const int g = idxTarget[d];
                                idxTarget[d] = (s == 0) ? -1 - g : m_res[d] + g;
                                idxSource = idxTarget;
                                if(bc[s][d] == TRANSMISSIVE)
                                {
                                    idxSource[d] = (s == 0) ? 0 : m_res[d] - 1;
                                }
                                else if(bc[s][d] == REFLECTIVE)
                                {
                                    idxSource[d] = (s == 0) ? g : m_res[d] - 1 - g;
                                }
                                else if(bc[s][d] == PERIODIC)
                                {
                                    idxSource[d] =  (s == 0) ? m_res[d] - 1 - g : g;
                                }
                                (*this)(idxTarget) = (*this)(idxSource);
                                if(bc[s][d] == REFLECTIVE)
                                {
                                    for(size_t v = 0; v < vecIdx.size(); ++v)
                                    {
                                        (*this)(idxTarget)[vecIdx[v][d]] *= -1.0;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        void save(const std::string& filename, const int nGhost=0)
        {
            assert(nGhost <= m_nGhost);
            const int nCellsInclGhost = GRIDDIM_TERM((m_res[0] + 2 * nGhost), * (m_res[1] + 2 * nGhost), * (m_res[2] + 2 * nGhost));
            m_buffer.resize(nCellsInclGhost * NVARS);
            size_t idx = 0;
            for(int i = -nGhost; i < m_res[0] + nGhost; ++i)
            {
                #if GRIDDIM >= 2
                for(int j = -nGhost; j < m_res[1] + nGhost; ++j)
                #endif
                {
                    #if GRIDDIM == 3
                    for(int k = -nGhost; k < m_res[2] + nGhost; ++k)
                    #endif
                    {
                        const std::array<REAL, NVARS>& data = (*this)(GRIDDIM_DECL(i, j, k));
                        std::copy_n(data.begin(), NVARS, &m_buffer[idx]);
                        idx += NVARS;
                    }
                }
            }
            std::ofstream out(filename, std::ios::binary);
            out.write(reinterpret_cast<const char*>(m_buffer.data()), m_buffer.size() * sizeof(REAL));
            out.close();
        }

        void load(const std::string& filename, const int nGhost=0)
        {
            assert(nGhost <= m_nGhost);
            const int nCellsInclGhost = GRIDDIM_TERM((m_res[0] + 2 * nGhost), * (m_res[1] + 2 * nGhost), * (m_res[2] + 2 * nGhost));
            m_buffer.resize(nCellsInclGhost * NVARS);
            std::ifstream in(filename, std::ios::binary);
            in.read(reinterpret_cast<char*>(m_buffer.data()), m_buffer.size() * sizeof(REAL));
            in.close();
            size_t idx = 0;
            for(int i = -nGhost; i < m_res[0] + nGhost; ++i)
            {
                #if GRIDDIM >= 2
                for(int j = -nGhost; j < m_res[1] + nGhost; ++j)
                #endif
                {
                    #if GRIDDIM == 3
                    for(int k = -nGhost; k < m_res[2] + nGhost; ++k)
                    #endif
                    {
                        std::array<REAL, NVARS>& data = (*this)(GRIDDIM_DECL(i, j, k));
                        std::copy_n(&m_buffer[idx], NVARS, data.begin());
                        idx += NVARS;
                    }
                }
            }
        }

    private:
        const std::array<int, GRIDDIM> m_res;
        const int m_nGhost;

        #if GRIDDIM == 1
            std::vector<std::array<REAL, NVARS>> m_data;
        #elif GRIDDIM == 2
            std::vector<std::vector<std::array<REAL, NVARS>>> m_data;
        #else  // GRIDDIM == 3
            std::vector<std::vector<std::vector<std::array<REAL, NVARS>>>> m_data;
        #endif

        std::vector<REAL> m_buffer;
};

template <int NVARS>
class Mesh : public DataArray<NVARS>
{
    public:
        #ifdef USE_RIGID
            static constexpr int SDFNGHOST = 1;
        #endif

        Mesh(const Geometry& geom, const int nGhost=0) : DataArray<NVARS>(geom.getRes(), nGhost), m_geom(geom) 
        #ifdef USE_RIGID
            , m_sdf(geom.getRes(), SDFNGHOST)
        #endif
        {
            
        }

        const Geometry& getGeometry() const 
        {
            return m_geom;
        }

        static Mesh<NVARS> createFromFile(const std::string& headerFilename, int& step, REAL& time, const int nGhost=0)
        {
            std::ifstream file(headerFilename);
            assert(file.is_open());
            std::string stepLine, timeLine, loLine, hiLine, resLine, NVARSLine, dataFilename;
            std::getline(file, stepLine);
            std::getline(file, timeLine);
            std::getline(file, loLine);
            std::getline(file, hiLine);
            std::getline(file, resLine);
            std::getline(file, NVARSLine);
            std::getline(file, dataFilename);
            file.close();
            std::istringstream stepISS(stepLine);
            std::istringstream timeISS(timeLine);
            std::istringstream loISS(loLine);
            std::istringstream hiISS(hiLine);
            std::istringstream resISS(resLine);
            std::istringstream NVARSISS(NVARSLine);
            stepISS >> step;
            timeISS >> time;
            std::array<REAL, GRIDDIM> lo, hi;
            std::array<int, GRIDDIM> res;
            for(int d = 0; d < GRIDDIM; ++d)
            {
                loISS >> lo[d];
                hiISS >> hi[d];
                resISS >> res[d];
            }
            int NVARSFromFile;
            NVARSISS >> NVARSFromFile;
            assert(NVARSFromFile == NVARS);
            const Geometry geom(lo, hi, res);
            Mesh<NVARS> mesh(geom, nGhost);
            mesh.load(dataFilename);
            return mesh;
        }

        void writeToFile(const std::string& headerFilename, const std::string& dataFilename, const int step, const REAL time)
        {
            std::ofstream file(headerFilename);
            assert(file.is_open());
            const std::array<REAL, GRIDDIM>& lo = m_geom.getLo();
            const std::array<REAL, GRIDDIM>& hi = m_geom.getHi();
            const std::array<int, GRIDDIM>& res = DataArray<NVARS>::getRes();
            file << step << std::endl;
            file << time << std::endl;
            file << GRIDDIM_TERM(lo[0], << " " << lo[1], << " " << lo[2]) << std::endl;
            file << GRIDDIM_TERM(hi[0], << " " << hi[1], << " " << hi[2]) << std::endl;
            file << GRIDDIM_TERM(res[0], << " " << res[1], << " " << res[2]) << std::endl;
            file << NVARS << std::endl;
            file << dataFilename << std::endl;
            file.close();
            this->save(dataFilename);
        }

        #ifdef USE_RIGID
            void setSDF(const REAL sdf, const std::array<int, GRIDDIM>& idx)
            {
                m_sdf(idx)[0] = sdf;
            }

            bool isRigid(GRIDDIM_DECL(const int i, const int j, const int k)) const
            {
                return m_sdf(GRIDDIM_DECL(i, j, k))[0] < 0;
            }

            bool isRigid(const std::array<int, GRIDDIM> idx) const
            {
                return m_sdf(idx)[0] < 0;
            }

            void fillInternalGhost(const std::vector<std::array<int, GRIDDIM>>& vecIdx,
                                   const REAL multiplier=1.1)
            {
                const std::array<int, GRIDDIM>& res = DataArray<NVARS>::getRes();
                const int nGhost = DataArray<NVARS>::getNGhost();
                const std::array<REAL, GRIDDIM>& dx = m_geom.getDx();
                const REAL cellSize = std::sqrt(GRIDDIM_TERM(dx[0] * dx[0],
                                                             + dx[1] * dx[1],
                                                             + dx[2] * dx[2]));
                #ifdef USE_OMP
                #pragma omp parallel for default(none) shared(vecIdx, multiplier, res, nGhost, dx, cellSize) schedule(dynamic)
                #endif
                for(int i = 0; i < res[0]; ++i)
                {
                    #if GRIDDIM >= 2
                    for(int j = 0; j < res[1]; ++j)
                    #endif
                    {
                        #if GRIDDIM == 3
                        for(int k = 0; k < res[1]; ++k)
                        #endif
                        {
                            const std::array<int, GRIDDIM> idx = {GRIDDIM_DECL(i, j, k)};
                            const REAL sdf = m_sdf(GRIDDIM_DECL(i, j, k))[0];
                            if(isRigid(idx) && sdf > -multiplier * nGhost * cellSize)
                            {
                                std::array<REAL, GRIDDIM> normal, pos, probePos;
                                std::array<int, GRIDDIM> offset;
                                m_geom.getPos(pos, idx);
                                for(int d = 0; d < GRIDDIM; ++d)
                                {
                                    offset = {GRIDDIM_DECL(0, 0, 0)};
                                    offset[d] = 1;
                                    const REAL sdfLo = m_sdf(GRIDDIM_DECL(i - offset[0],
                                                                          j - offset[1],
                                                                          k - offset[2]))[0];
                                    const REAL sdfHi = m_sdf(GRIDDIM_DECL(i + offset[0],
                                                                          j + offset[1],
                                                                          k + offset[2]))[0];
                                    normal[d] = 0.5 * (sdfHi - sdfLo) / dx[d];
                                    probePos[d] = pos[d] - sdf * normal[d] * multiplier;
                                }
                                std::array<REAL, GRIDDIM> probeIdx;
                                m_geom.getIdx(probeIdx, probePos);
                                std::array<REAL, NVARS> probeData = (*this)(probeIdx);
                                std::array<REAL*, GRIDDIM> vec;
                                for(size_t v = 0; v < vecIdx.size(); ++v)
                                {
                                    vec = {GRIDDIM_DECL(&probeData[vecIdx[v][0]], 
                                                        &probeData[vecIdx[v][1]], 
                                                        &probeData[vecIdx[v][2]])};
                                    const REAL dot = GRIDDIM_TERM(normal[0] * (*vec[0]), 
                                                                  + normal[1] * (*vec[1]), 
                                                                  + normal[2] * (*vec[2]));
                                    for(int d = 0; d < GRIDDIM; ++d)
                                    {
                                        (*vec[d]) -= 2 * normal[d] * dot;
                                    }
                                }
                                for(int v = 0; v < NVARS; ++v)
                                {
                                    (*this)(GRIDDIM_DECL(i, j, k))[v] = probeData[v];
                                }
                            }
                        }
                    }
                }
            }

            void readSDFFromFile(const std::string& headerFilename)
            {
                const std::array<REAL, GRIDDIM>& lo = m_geom.getLo();
                const std::array<REAL, GRIDDIM>& hi = m_geom.getHi();
                const std::array<int, GRIDDIM>& res = DataArray<NVARS>::getRes();
                std::ifstream file(headerFilename);
                assert(file.is_open());
                std::string loLine, hiLine, resLine, dataFilename;
                std::getline(file, loLine);
                std::getline(file, hiLine);
                std::getline(file, resLine);
                std::getline(file, dataFilename);
                file.close();
                std::istringstream loISS(loLine);
                std::istringstream hiISS(hiLine);
                std::istringstream resISS(resLine);
                REAL loFromFile, hiFromFile;
                int resFromFile;
                for(int d = 0; d < GRIDDIM; ++d)
                {
                    loISS >> loFromFile;
                    hiISS >> hiFromFile;
                    resISS >> resFromFile;
                    assert(std::fabs(loFromFile - lo[d]) < 1e-16);
                    assert(std::fabs(hiFromFile - hi[d]) < 1e-16);
                    assert(resFromFile == res[d]);
                }
                m_sdf.load(dataFilename, 1);
            }

            void writeSDFToFile(const std::string& headerFilename, const std::string& dataFilename)
            {
                std::ofstream file(headerFilename);
                assert(file.is_open());
                const std::array<REAL, GRIDDIM>& lo = m_geom.getLo();
                const std::array<REAL, GRIDDIM>& hi = m_geom.getHi();
                const std::array<int, GRIDDIM>& res = DataArray<NVARS>::getRes();
                file << GRIDDIM_TERM(lo[0], << " " << lo[1], << " " << lo[2]) << std::endl;
                file << GRIDDIM_TERM(hi[0], << " " << hi[1], << " " << hi[2]) << std::endl;
                file << GRIDDIM_TERM(res[0], << " " << res[1], << " " << res[2]) << std::endl;
                file << dataFilename << std::endl;
                file.close();
                m_sdf.save(dataFilename, 1);
            }
        #endif

    private:
        const Geometry m_geom;

        #if GRIDDIM == 1
            std::vector<std::array<REAL, NVARS>> m_data;
        #elif GRIDDIM == 2
            std::vector<std::vector<std::array<REAL, NVARS>>> m_data;
        #else  // GRIDDIM == 3
            std::vector<std::vector<std::vector<std::array<REAL, NVARS>>>> m_data;
        #endif

        #ifdef USE_RIGID
            DataArray<1> m_sdf;
        #endif
};

#endif