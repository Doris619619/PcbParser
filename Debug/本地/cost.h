#ifndef COST_H
#define COST_H

#include "grid.h"
#include <vector>
#include <string>
#include "Separate3.h"


class CostUpdater {
public:
    CostUpdater(Grid& grid);
    void updateCost(int layerId,int start_x, int start_y, int end_x, int end_y);
    void outputWeightsToFile(int layerId,int start_x, int start_y, int end_x, int end_y, const std::string& outputFile);
private:
    Grid& grid;  
    int directions[4][2] = { {-1, 0}, {1, 0}, {0, -1}, {0, 1} };  

};


#endif
#pragma once
