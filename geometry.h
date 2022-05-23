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
    double function(double x);
};

class segment : public graph
{
public:
    segment(Point A, Point B);
    double getLength();
    double getAngle(segment* line)
    {
        Point A = line->getVertices().first;
        Point B = line->getVertices().second;
        std::pair<double,double> firstV;
        std::pair<double,double> secondV;
        if(m_A.x == A.x && m_A.y == A.y)
        {
            firstV = std::make_pair(m_B.x-m_A.x,m_B.y-m_A.y);
            secondV = std::make_pair(B.x-A.x,B.y-A.y);
        }
        else if(m_B.x == B.x && m_B.y == B.y)
        {
            firstV = std::make_pair(m_A.x-m_B.x,m_A.y-m_B.y);
            secondV = std::make_pair(A.x-B.x,A.y-B.y);
        }
        else if(m_A.x == B.x && m_A.y == B.y)
        {
            firstV = std::make_pair(m_B.x-m_A.x,m_B.y-m_A.y);
            secondV = std::make_pair(A.x-B.x,A.y-B.y);
        }
        else if(m_B.x == A.x && m_B.y == A.y)
        {
            firstV = std::make_pair(m_A.x-m_B.x,m_A.y-m_B.y);
            secondV = std::make_pair(B.x-A.x,B.y-A.y);
        }
        else{
            std::cout<<"does not intersect"<<std::endl;
            return 0;
        }

        double tmp = firstV.first*secondV.first + firstV.second*secondV.second;
        return acos(tmp/(m_length*line->getLength()));
    }
    std::pair<Point,Point> getVertices(){return std::make_pair(m_A,m_B);}
private:
    double m_length;
    Point m_A; Point m_B;
};

class figure{
public:
    void setFrequency(double frequency);

    std::vector<Point> getGraph();
    Point getCenter();

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

    double getSquare();
};
class quadrilateral : public figure{
public:
    quadrilateral(Point A,Point B, Point C, Point D);
    double getSquare();
};
#endif //GEOMETRY_LIBRARY_H
