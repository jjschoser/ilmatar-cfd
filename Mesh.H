#ifndef MESH_H
#define MESH_H

#include <array>
#include <cassert>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Macros.H"

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

        std::array<REAL, NVARS>& operator()(const std::array<int, GRIDDIM>& idx)
        {
            return (*this)(GRIDDIM_DECL(idx[0], idx[1], idx[2]));
        }

        const std::array<REAL, NVARS>& operator()(const std::array<int, GRIDDIM>& idx) const
        {
            return (*this)(GRIDDIM_DECL(idx[0], idx[1], idx[2]));
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
};

template <int NVARS>
class Mesh : public DataArray<NVARS>
{
    public:
        Mesh(const Geometry& geom, const int nGhost=0) : DataArray<NVARS>(geom.getRes(), nGhost), m_geom(geom) 
        {
            
        }

        const Geometry& getGeometry() const 
        {
            return m_geom;
        }

        static Mesh<NVARS> createFromFile(const std::string& filename, int& step, REAL& time, const int nGhost=0)
        {
            std::ifstream file(filename);
            assert(file.is_open());
            std::string stepLine, timeLine, loLine, hiLine, resLine, NVARSLine;
            std::getline(file, stepLine);
            std::getline(file, timeLine);
            std::getline(file, loLine);
            std::getline(file, hiLine);
            std::getline(file, resLine);
            std::getline(file, NVARSLine);
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
            std::string dataLine;
            int counter = 0;
            while(std::getline(file, dataLine))
            {
                const int i = counter % res[0];
                #if GRIDDIM >= 2
                const int j = (counter / res[0]) % res[1];
                #endif
                #if GRIDDIM == 3
                const int k = counter / (res[0] * res[1]);
                #endif
                std::istringstream dataISS(dataLine);
                for(int v = 0; v < NVARS; ++v)
                {
                    dataISS >> mesh({GRIDDIM_DECL(i, j, k)})[v];
                }
                ++counter;
            }
            return mesh;
        }

        void writeToFile(const std::string& filename, const int step, const REAL time)
        {
            std::ofstream file(filename);
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
            #if GRIDDIM == 3
            for(int k = 0; k < res[2]; ++k)
            #endif
            {
                #if GRIDDIM >= 2
                for(int j = 0; j < res[1]; ++j)
                #endif
                {
                    for(int i = 0; i < res[0]; ++i)
                    {
                        const std::array<REAL, NVARS> data = (*this)({GRIDDIM_DECL(i, j, k)});
                        for(int v = 0; v < NVARS; ++v)
                        {
                            file << data[v] << " ";
                        }
                        file << std::endl;
                    }
                }
            }
            file.close();
        }

    private:
        const Geometry m_geom;

        #if GRIDDIM == 1
            std::vector<std::array<REAL, NVARS>> m_data;
        #elif GRIDDIM == 2
            std::vector<std::vector<std::array<REAL, NVARS>>> m_data;
        #else  // GRIDDIM == 3
            std::vector<std::vector<std::vector<std::array<REAL, NVARS>>>> m_data;
        #endif
};

#endif