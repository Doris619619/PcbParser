#pragma once
#ifndef CALCULATOR_H
#define CALCULATOR_H

#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Footprint�ṹ�壨�����ڲ�������
struct Pad {
    std::string name;
    double x;
    double y;
};

struct Footprint {
    std::string name;
    double originX;
    double originY;
    double angle;
    std::vector<Pad> pads;
    Footprint() : originX(0), originY(0), angle(0) {}
};

// �µĽṹ�壬�洢footprint��pad�ľ���������Ϣ�����ڷ��ؽ����
struct FootprintPadAbsolute {
    std::string footprintName;
    std::string padName;
    double footprintOriginX;
    double footprintOriginY;
    double footprintAngle;
    double padAbsoluteX;
    double padAbsoluteY;
};

std::vector<FootprintPadAbsolute> calculateFootprintCoordinates(const std::string& inputFile);

#endif // CALCULATOR_H