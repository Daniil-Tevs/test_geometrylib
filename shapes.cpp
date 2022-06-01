#include "geometry.h"

double figure::getPerimeter(){
    double sum = 0;
    for(auto& i : m_sides)
        sum += i->getLength();
    return sum;
}

Point figure::getCenter(){return m_center;}
std::string figure::getType(){return m_type;}
std::vector<Point> figure::getPoints(){return m_points;}
std::vector<Point> figure::getVertices(){return m_vertices;}
std::vector<segment*> figure::getSides() const {return m_sides;}

triangle::triangle(Point A, Point B, Point C)
{
    m_type = "triangle";
    m_vertices.push_back(A); m_vertices.push_back(B); m_vertices.push_back(C);
    m_center = { (A.x+B.x+C.x)/3,(A.y+B.y+C.y)/3};
    m_sides.push_back(new segment(A,B));
    m_sides.push_back(new segment(B,C));
    m_sides.push_back(new segment(A,C));
    for(auto& i : m_sides)
        for(auto&j : i->getGraph())
            m_points.push_back(j);
    for(auto& i : m_vertices)
        m_points.push_back(i);
    std::sort(m_points.begin(), m_points.end());
    for(int i=0;i<m_points.size();i++)
        for(int j=i+1;j<m_points.size();j++) {
            if (m_points[i] == m_points[j]) {
                m_points.erase(m_points.begin() + j);
            }
        }
    for(auto& i : m_points)
        m_base.push_back(i);
}
double triangle::getSquare()
{
    double AB = m_sides[0]->getLength(), BC = m_sides[1]->getLength(),AC = m_sides[2]->getLength();
    double p = (AB+BC+AC)/2;
    return sqrt(p*(p-AB) * (p-BC) *(p-AC));
}

quadrilateral::quadrilateral(Point A,Point B, Point C, Point D)
{
    m_type = "quadrilateral";
    m_vertices.push_back(A); m_vertices.push_back(B);
    m_vertices.push_back(C); m_vertices.push_back(D);
    double tmp1=abs((B.x-A.x)*(C.y-A.y)-(C.x-A.y)*(B.y-A.y))/2;
    double tmp2=abs((C.x-A.x)*(D.y-A.y)-(D.x-A.y)*(C.y-A.y))/2;
    Point tmpX1={(A.x+B.x+C.x)/3,(A.y+B.y+C.y)/3};
    Point tmpX2= {(A.x + C.x + D.x) / 3,(A.y+C.y+D.y)/3};
    m_center = { (tmp1*tmpX1.x+tmp2*tmpX2.x)/(tmp1+tmp2),(tmp1*tmpX1.y+tmp2*tmpX2.y)/(tmp1+tmp2)};

    m_sides.push_back(new segment(A,B));
    m_sides.push_back(new segment(B,C));
    m_sides.push_back(new segment(C,D));
    m_sides.push_back(new segment(D,A));
    for(auto& i : m_vertices)
        m_points.push_back(i);
    for(auto& i : m_sides)
        for(auto& j : i->getGraph())
            m_points.push_back(j);
    std::sort(m_points.begin(), m_points.end());
    for(int i=0;i<m_points.size();i++)
        for(int j=i+1;j<m_points.size();j++) {
            if (m_points[i] == m_points[j]) {
                m_points.erase(m_points.begin() + j);
            }
        }
    for(auto& i : m_points)
        m_base.push_back(i);
}
double quadrilateral::getSquare(){
    double AB = m_sides[0]->getLength(), BC = m_sides[1]->getLength(),CD = m_sides[2]->getLength(),AD = m_sides[3]->getLength();
    double p = (AB+BC+CD+AD)/2;
    double cosAngle = cos((m_sides[0]->getAngle(m_sides[1])+m_sides[2]->getAngle(m_sides[3]))/2);
    return sqrt((p-AB) * (p-BC) *(p-CD)*(p-AD)-AB * BC*CD*AD* pow(cosAngle,2));
}

circle::circle(Point center,double R){
    m_type = "circle";
    m_center = center;
    m_radius = R;
    double tmp;
    for(double i = m_center.x - R;i<=m_center.x+R;i+=m_frequency)
    {
        tmp = sqrt(m_radius*m_radius - pow((i-m_center.x),2));
        m_points.push_back({i,m_center.y - tmp});
        m_points.push_back({i,m_center.y+tmp});
    }
    m_points.push_back({m_center.x+R * cos(0),m_center.y+R * sin(0)});
    m_points.push_back({m_center.x+R * cos(acos(-1)/2),m_center.y+R * sin(acos(-1)/2)});
    m_points.push_back({m_center.x+R * cos(acos(-1)),m_center.y+R * sin(acos(-1))});
    m_points.push_back({m_center.x+R * cos(3*acos(-1)/2),m_center.y+R * sin(3*acos(-1)/2)});
//    double tmpX,tmpY;
//    for(double fi = 0;fi<=2*acos(-1);fi+=m_frequency)
//    {
//        tmpX = R*cos(fi);
//        tmpY = R*sin(fi);
//        m_points.push_back({tmpX, tmpY});
//    }

    for(int i=0;i<m_points.size();i++)
        for(int j=i+1;j<m_points.size();j++) {
            if (m_points[i].x == m_points[j].x && m_points[i].y == m_points[j].y) {
                m_points.erase(m_points.begin() + j);
                i--;break;
            }
        }
    for(auto& i : m_points)
        m_base.push_back(i);
}
double circle::getSquare(){return acos(-1)*m_radius*m_radius;}
double circle::getPerimeter(){
    return 2*acos(-1)*m_radius;
}
double circle::getRadius() const{return m_radius;}
std::pair<double,double> circle::function(double x)
{
    double tmp = sqrt(m_radius*m_radius - pow((x-m_center.x),2));
    return std::make_pair(m_center.y+tmp,m_center.y - tmp);
}

