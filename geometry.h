#ifndef GEOMETRY_LIBRARY_H
#define GEOMETRY_LIBRARY_H
#include <iostream>
#include <cmath>
#include <vector>
struct Point{
    double x;
    double y;
};
class graph
{
public:
    graph( Point A, Point B, std::pair<double,double> D_x);
    ~graph() = default;
    std::vector<Point> getGraph();
    void setFrequency(double frequency);
protected:
    Point M_0;
    std::vector<Point> m_points;
    double m_frequency;
    std::pair<double,double> collinear;
    std::pair<double,double> m_D_x;
    double function(double x) const;
};

class segment : public graph
{
public:
    segment(Point A, Point B);
    double getLength() const;
    double getAngle(segment* line) const;
    std::pair<Point,Point> getVertices();
private:
    double m_length;
    Point m_A; Point m_B;
};

class figure{
public:
    void setFrequency(double frequency);

    std::vector<Point> getGraph();
    Point getCenter();
    std::vector<Point> getVertices();

    virtual double getSquare() = 0;
protected:
    figure() = default;
    ~figure() = default;
    std::vector<Point> m_points;
    std::vector<segment*> m_sides;
    std::vector<Point> m_vertices;
    Point m_center;
    double m_frequency=1;
};

class triangle : public figure{
public:
    triangle(Point A, Point B, Point C);
    virtual ~triangle() = default;

    double getSquare() override;
};

class quadrilateral : public figure{
public:
    quadrilateral(Point A,Point B, Point C, Point D);
    double getSquare() override;
};
#endif //GEOMETRY_LIBRARY_H
