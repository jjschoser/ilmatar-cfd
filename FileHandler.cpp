#include "FileHandler.H"

#include <filesystem>

// Contents of this file were written with the help of ChatGPT

std::string getDataFilename(const std::string& headerFilename)
{
    std::filesystem::path path(headerFilename);
    path.replace_filename(path.stem().string() + "Data.dat");
    return path.string();
}

std::string removePath(const std::string& filename)
{
    std::filesystem::path path(filename);
    return path.filename().string();
}

std::string addPath(const std::string& srcFilename, const std::string& dstFilename)
{
    std::filesystem::path srcPath(srcFilename);
    std::filesystem::path newPath = srcPath.parent_path() / dstFilename;
    return newPath.string();
}

std::string addStepCounter(const std::string& filename, const int step)
{
    std::filesystem::path path(filename);
    std::string dir = path.parent_path().string();
    std::string stem = path.stem().string();
    std::string extension = path.extension().string();
    return dir + "/" + stem + std::to_string(step) + extension;
}

void createDirForFile(const std::string& filename)
{
    std::filesystem::path path(filename);
    std::filesystem::path dir = path.parent_path();
    if(!std::filesystem::exists(dir))
    {
        std::filesystem::create_directory(dir);
    }
}
