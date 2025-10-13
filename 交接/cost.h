#ifndef COST_H
#define COST_H

#include "grid.h"
#include <vector>
#include <string>

class CostUpdater {
public:
    CostUpdater(Grid& grid);
    void updateCost();
    void outputWeightsToFile(const std::string& outputFile);
private:
    Grid& grid;  // 引用传递Grid对象
    int directions[4][2] = { {-1, 0}, {1, 0}, {0, -1}, {0, 1} };  // 上下左右四个方向

};

#endif
#pragma once
