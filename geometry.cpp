#include "geometry.h"

graph::graph( Point A, Point B, std::pair<double,double> D_x)
{
    m_D_x = D_x;
    m_frequency = 1; M_0 = {0,0};
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
        for(double i = D_x.first ;i<=D_x.second;i+=m_frequency)
            m_points.push_back({double(i),function(i)});
    else
        for(double i = D_x.first ;i<=D_x.second;i+=m_frequency)
            m_points.push_back({0,i});
}

std::vector<Point> graph::getGraph(){return m_points;}

void graph::setFrequency(double frequency){
    m_frequency = abs(frequency);
    m_points.clear();
    if(collinear.first!=0)
        for(double i = m_D_x.first ;i<=m_D_x.second;i+=m_frequency)
            m_points.push_back({double(i),function(i)});
    else
        for(double i = m_D_x.first ;i<=m_D_x.second;i+=m_frequency)
            m_points.push_back({0,i});
}

double graph::function(double x) const {return collinear.second * (x-M_0.x)/collinear.first + M_0.y;}

segment::segment(Point A, Point B) : graph(A,B,{(A.x<=B.x)?A.x:B.x,(B.x<=A.x)?A.x:B.x}) {
    m_length = sqrt(pow(A.x - B.x,2)+pow(A.y - B.y,2));
    m_A = A; m_B = B;
}
double segment::getLength() const {
    return m_length;}

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


std::vector<Point> figure::getGraph(){return m_points;}
void figure::setFrequency(double frequency)
{
    m_frequency = abs(frequency);
    for(auto& i : m_sides)
        i->setFrequency(m_frequency);
    m_points.clear();
    for(auto& i : m_sides)
        for(auto&j : i->getGraph())
            m_points.push_back(j);
    for(int i=0;i<m_points.size();i++)
        for(int j=i+1;j<m_points.size();j++) {
            if (m_points[i].x >= m_points[j].x)
                std::swap(m_points[i], m_points[j]);
            if (m_points[i].x == m_points[j].x && m_points[i].y == m_points[j].y) {
                m_points.erase(m_points.begin() + j);
                i--;break;
            }
        }
}
Point figure::getCenter(){return m_center;}
std::vector<Point> figure::getVertices(){return m_vertices;}

triangle::triangle(Point A, Point B, Point C)
{
    m_vertices.push_back(A);
    m_vertices.push_back(B);
    m_vertices.push_back(C);
    m_center = { 1/3*(A.x+B.x+C.x),1/3*(A.y+B.y+C.y)};
    m_sides.push_back(new segment(A,B));
    m_sides.push_back(new segment(B,C));
    m_sides.push_back(new segment(A,C));
    for(auto& i : m_sides)
        for(auto&j : i->getGraph())
            m_points.push_back(j);
    for(int i=0;i<m_points.size();i++)
        for(int j=i+1;j<m_points.size();j++) {
            if (m_points[i].x >= m_points[j].x)
                std::swap(m_points[i], m_points[j]);
            if (m_points[i].x == m_points[j].x && m_points[i].y == m_points[j].y) {
                m_points.erase(m_points.begin() + j);
                i--;break;
            }
        }
}
double triangle::getSquare()
{
    double AB = m_sides[0]->getLength(), BC = m_sides[1]->getLength(),AC = m_sides[2]->getLength();
    double p = (AB+BC+AC)/2;
    return sqrt(p*(p-AB) * (p-BC) *(p-AC));
}

quadrilateral::quadrilateral(Point A,Point B, Point C, Point D)
{
    m_vertices.push_back(A);
    m_vertices.push_back(B);
    m_vertices.push_back(C);
    m_vertices.push_back(D);
    double tmp1=1/2*abs((B.x-A.x)*(C.y-A.y)-(C.x-A.y)*(B.y-A.y));
    double tmp2=1/2*abs((C.x-A.x)*(D.y-A.y)-(D.x-A.y)*(C.y-A.y));
    Point tmpX1={(A.x+B.x+C.x)/3,(A.y+B.y+C.y)/3};
    Point tmpX2= {(A.x + C.x + D.x) / 3,(A.y+C.y+D.y)/3};
    m_center = { (tmp1*tmpX1.x+tmp2*tmpX2.x)/(tmp1+tmp2),(tmp1*tmpX1.y+tmp2*tmpX2.y)/(tmp1+tmp2)};

    m_sides.push_back(new segment(A,B));
    m_sides.push_back(new segment(B,C));
    m_sides.push_back(new segment(C,D));
    m_sides.push_back(new segment(A,D));
    for(auto& i : m_sides)
        for(auto&j : i->getGraph())
            m_points.push_back(j);
    for(int i=0;i<m_points.size();i++)
        for(int j=i+1;j<m_points.size();j++) {
            if (m_points[i].x >= m_points[j].x)
                std::swap(m_points[i], m_points[j]);
            if (m_points[i].x == m_points[j].x && m_points[i].y == m_points[j].y) {
                m_points.erase(m_points.begin() + j);
                i--;break;
            }
        }
}
double quadrilateral::getSquare(){
    double AB = m_sides[0]->getLength(), BC = m_sides[1]->getLength(),CD = m_sides[2]->getLength(),AD = m_sides[3]->getLength();
    double p = (AB+BC+CD+AD)/2;
    double cosAngle = cos((m_sides[0]->getAngle(m_sides[1])+m_sides[2]->getAngle(m_sides[3]))/2);
    return sqrt((p-AB) * (p-BC) *(p-CD)*(p-AD)-AB * BC*CD*AD* pow(cosAngle,2));
}