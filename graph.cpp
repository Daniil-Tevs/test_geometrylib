#include "geometry.h"

graph::graph( Point A, Point B, std::pair<double,double> D_x)
{
    m_D_x = D_x;
    M_0 = {0,0};
    if(A.x==B.x && A.y==B.y) { std::cerr << "Error: equal points" << std::endl; collinear = {0,0};}
    if(B.x>A.x) {
        M_0 = A;
        collinear = std::make_pair(B.x - A.x,B.y-A.y);
    }
    else
    {
        M_0 = B;
        collinear = std::make_pair(A.x - B.x,A.y-B.y);
    }
    if(collinear.first!=0)
        for(double i = D_x.first ;i<D_x.second;i+=m_frequency)
            m_points.push_back({double(i),function(i)});
    else
    {
        m_D_x.first = (A.y<=B.y)?A.y:B.y;
        m_D_x.second = (A.y>=B.y)?A.y:B.y;
        for(double i = m_D_x.first ;i<m_D_x.second;i+=m_frequency)
            m_points.push_back({A.x,i});
    }
    for(auto& i : m_points)
        m_base.push_back(i);
}

std::vector<Point> graph::getGraph(){return m_points;}
std::vector<Point> graph::getBasePoints() {return m_base;}
double segment::getLength() const {return m_length;}

void graph::setFrequency(double frequency){
    m_frequency = abs(frequency);
    m_points.clear();
    if(collinear.first!=0)
        for(double i = m_D_x.first ;i<=m_D_x.second;i+=m_frequency)
            m_points.push_back({double(i),function(i)});
    else
        for(double i = m_D_x.first ;i<=m_D_x.second;i+=m_frequency)
            m_points.push_back({M_0.x,i});
}

double graph::function(double x) const {return collinear.second * (x-M_0.x)/collinear.first + M_0.y;}

segment::segment(Point A, Point B) : graph(A,B,{(A.x<=B.x)?A.x:B.x,(B.x<=A.x)?A.x:B.x}) {
    m_length = sqrt(pow(A.x - B.x,2)+pow(A.y - B.y,2));
    m_A = A; m_B = B;
}

double segment::getAngle(segment* line) const
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

std::pair<Point,Point> segment::getVertices(){return std::make_pair(m_A,m_B);}

std::vector<Point> segment::getIntersectionPoints(const segment& segment) const
{
    std::vector<Point> inPoints;
    for (auto &i: m_points)
        for (auto &j: segment.m_points)
            if (abs(i.x - j.x) < m_accuracy && abs(i.y - j.y) < m_accuracy)
                inPoints.push_back(i);
    std::reverse(inPoints.begin(), inPoints.end());
    if(!inPoints.empty()) {
        int t=0;
        while(t<2) {
            for (int i = 0; i < inPoints.size() - 1; i++)
                for (int j = i + 1; j < inPoints.size(); j++) {
                    if (abs(inPoints[i].x - inPoints[j].x) < m_accuracy &&
                        abs(inPoints[i].y - inPoints[j].y) < m_accuracy) {
                        inPoints.erase(inPoints.begin() + j);
                        break;
                    }
                    if (round(inPoints[i].x) == round(inPoints[j].x) &&
                        round(inPoints[i].y) == round(inPoints[j].y)) {
                        int tmpX = round(inPoints[i].x);
                        int tmpY = round(inPoints[i].y);
                        if(abs(tmpX-inPoints[i].x)<=m_accuracy*5)
                            inPoints[i].x = round(inPoints[i].x);
                        if(abs(tmpY -  inPoints[i].y)<=m_accuracy*5)
                            inPoints[i].y = round(inPoints[i].y);
                        inPoints.erase(inPoints.begin() + j);
                        break;
                    }
                }
            t++;
        }
    }
    return inPoints;
}
