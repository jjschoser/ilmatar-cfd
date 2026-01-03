#ifndef FILE_HANDLER_H
#define FLUX_HANDLER_H

#include <string>

#include "Macros.H"

std::string getDataFilename(const std::string& headerFilename);
std::string removePath(const std::string& filename);
std::string addPath(const std::string& srcFilename, const std::string& dstFilename);
std::string addStepCounter(const std::string& filename, const int step);
void createDirForFile(const std::string& filename);

#endif