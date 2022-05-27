#ifndef GEOMETRY_LIBRARY_H
#define GEOMETRY_LIBRARY_H
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#define basic_frequency 0.001
#define basic_accuracy 0.001

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

class graph
{
public:
    graph( Point A, Point B, std::pair<double,double> D_x);
    ~graph() = default;
    std::vector<Point> getGraph();
    void setFrequency(double frequency);
    std::vector<Point> getBasePoints();
protected:
    Point M_0;
    std::vector<Point> m_points;
    std::vector<Point> m_base;
    double m_frequency = basic_frequency;
    double m_accuracy =basic_frequency;
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
    std::vector<Point> getIntersectionPoints(const segment& segment) const;
private:
    double m_length;
    Point m_A; Point m_B;
};

/// CLass figure - is parent class for shapes
// setFrequency() - setting the frequency of points
// getGraph() - returns the points of the shape
// getCenter() - returns the center of this shape
// getSquare() - returns square of this shape
// getPerimeter() - returns perimeter of shape
// m_points - an array of shape points that the user can change
// m_base - array of points of the figure that are the main ones for the initial frequency - basic_frequency

class figure{
public:
    virtual void setFrequency(double frequency);
    void setAccuracy(double accuracy);
    std::vector<Point> getGraph();
    Point getCenter();
    std::vector<Point> getVertices();
    std::string getType();
    std::vector<Point> getBasePoint() const {return m_base;}
    std::vector<segment*> getSides() const {return m_sides;}

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

    void setFrequency(double frequency) override;
    double getSquare() override;
    double getPerimeter() override;
private:
    double m_radius;
};

// Можно оставить
double getLength(Point A,Point B);
double getAngle(segment* AB,segment* CD);
double getAngleSegments(Point A,Point B,Point C,Point D);
double getPerimeter(Point A,Point B,Point C);
double getPerimeter(Point A,Point B,Point C, Point D);
double getSquare(Point A,Point B,Point C);
double getSquare(Point A,Point B,Point C,Point D);
///До сюда

std::vector<Point> getIntersectionPoints(segment* first, segment* segment);
std::vector<Point> getIntersectionPoints(triangle* first, triangle* second);
double getIntersectionSquare(triangle* first, triangle* second);


#endif //GEOMETRY_LIBRARY_H
