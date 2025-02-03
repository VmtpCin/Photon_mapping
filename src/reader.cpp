#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "reader.h"

std::vector<Object> read(char *filename) {
    std::vector<Point3> pnts;
    std::vector<Object> objs;

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return objs;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string prefix;
        ss >> prefix;

        if (prefix == "v") {
            Point3 pnt;
            ss >> pnt[0] >> pnt[1] >> pnt[2];
            pnts.push_back(pnt);
        } else if (prefix == "vt") {
        } else if (prefix == "vn") {
        } else if (prefix == "f") {
            
        }
    }

    file.close();

    return objs;
}