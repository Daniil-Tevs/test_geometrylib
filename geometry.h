#ifndef GEOMETRY_LIBRARY_H
#define GEOMETRY_LIBRARY_H
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#define basic_frequency 0.001
#define basic_accuracy 0.002999

///Надо добавит проверку для треугольника принадл одной прямой
///Проверку для четырехугольника
///настроить проверку на пересечение
///добавить getIntr для круга и прямоугольников

struct Point{
    double x;
    double y;
    Point operator+(const Point& A);
    Point operator-(const Point& A);
    bool operator==(const Point& A) const;
    bool operator>(const Point& A) const;
    bool operator<(const Point& A) const;
    friend std::ostream& operator<<(std::ostream &os, const Point& A) {
        return os  << '(' << A.x << ';' << A.y << ')';
    }
};
double getLengthArc(double radius)
{return 2*acos(-1)*pow(radius,2);}
double getSquare(double radius)
{return acos(-1)*pow(radius,2);}
class graph
{
public:
    graph( Point A, Point B, std::pair<double,double> D_x);
    ~graph() = default;
    std::vector<Point> getGraph();
    double functionY(double x) const;
    double functionX(double y) const;
    std::pair<double,double> getCollinear();
protected:
    Point M_0;
    std::vector<Point> m_points;
    double m_frequency = basic_frequency;
    double m_accuracy =basic_frequency;
    std::pair<double,double> collinear;
    std::pair<double,double> m_D_x;

};

class segment : public graph
{
public:
    segment(Point A, Point B);
    double getLength() const;
    double getAngle(segment* line) const;
    std::pair<Point,Point> getVertices();
    std::vector<Point> getIntersectionPoints(const segment& segment) const;
private:
    double m_length;
    Point m_A; Point m_B;
};

/// CLass figure - is parent class for shapes
// getPoints() - returns the points of the shape
// getCenter() - returns the center of this shape
// getSquare() - returns square of this shape
// getPerimeter() - returns perimeter of shape
// m_points - an array of shape points that the user can change
// m_base - array of points of the figure that are the main ones for the initial frequency - basic_frequency

class figure{
public:
    std::vector<Point> getPoints();
    Point getCenter();
    std::vector<Point> getVertices();
    std::vector<segment*> getSides() const;
    std::string getType();
    virtual double getSquare() = 0;
    virtual double getPerimeter();
protected:
    figure() = default;
    ~figure() = default;
    std::vector<Point> m_points;
    std::vector<Point> m_base;
    std::vector<segment*> m_sides;
    std::vector<Point> m_vertices;
    Point m_center;
    double m_frequency=basic_frequency;
    double m_accuracy = basic_frequency;
    std::string m_type;
};

///Class triangle
class triangle : public figure{
public:
    triangle(Point A, Point B, Point C);
    virtual ~triangle() = default;

    double getSquare() override;
};
///Class quadrilateral
class quadrilateral : public figure{
public:
    quadrilateral(Point A,Point B, Point C, Point D);
    virtual ~quadrilateral() = default;
    double getSquare() override;
};

class circle : public figure{
public:
    circle(Point center,double R);
    virtual ~circle() = default;

    double getSquare() override;
    double getPerimeter() override;
    double getRadius() const;
    std::pair<double,double> function(double x);
private:
    double m_radius;
};

double getLength(Point A,Point B);
double getAngle(segment* AB,segment* CD);
double getAngleSegments(Point A,Point B,Point C,Point D);
double getPerimeter(Point A,Point B,Point C);
double getPerimeter(Point A,Point B,Point C, Point D);
double getSquare(Point A,Point B,Point C);
double getSquare(Point A,Point B,Point C,Point D);

std::vector<Point> getIntersectionPoints(segment* AB, segment* CD);
//Добавить совпадение фигуры и сторон для нахождения площади. Поправить нахождение площади 5-6 угол
std::vector<Point> getIntersectionPoints(triangle* first, triangle* second);
//Добавить совпадение фигуры и сторон для нахождения площади
//+случай 7 и 8 угольника
std::vector<Point> getIntersectionPoints(quadrilateral* first, quadrilateral* second);
std::vector<Point> getIntersectionPoints(circle* first, circle* second);
std::vector<Point> getIntersectionPoints(segment* AB, circle* shape);
std::vector<Point> getIntersectionPoints(triangle* first, circle* second);
std::vector<Point> getIntersectionPoints(quadrilateral* first, circle* second);

double getIntersectionSquare(triangle* first, triangle* second);
double getIntersectionSquare(quadrilateral* first, quadrilateral* second);
double getIntersectionSquare(circle* first, circle* second);
double getIntersectionSquare(triangle* first, quadrilateral* second);
double getIntersectionSquare(triangle* first, circle* second);
double getIntersectionSquare(quadrilateral* first, circle* second);
#endif //GEOMETRY_LIBRARY_H
