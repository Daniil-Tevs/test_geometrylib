#include "geometry.h"

Point Point::operator+(const Point& A)
{
    this->x += A.x; this->y += A.y;
    return *this;
}
Point Point::operator-(const Point& A){
    this->x -= A.x; this->y -= A.y;
    return *this;
}
bool Point::operator==(const Point& A) const {return (A.x == this->x) && (A.y == this->y);}
bool Point::operator>(const Point& A) const {return (this->x>A.x);}
bool Point::operator<(const Point& A) const {return (this->x < A.x);}

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
            m_points.push_back({M_0.x,i});
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


std::vector<Point> figure::getGraph(){return m_points;}
double figure::getPerimeter(){
    double sum = 0;
    for(auto& i : m_sides)
        sum += i->getLength();
    return sum;
}
void figure::setFrequency(double frequency)
{
    m_frequency = abs(frequency);
    for(auto& i : m_sides)
        i->setFrequency(m_frequency);
    m_points.clear();
    for(auto& i : m_sides)
        for(auto& j : i->getGraph())
            m_points.push_back(j);
    for(auto& i : m_vertices)
        m_points.push_back(i);
    for(int i=0;i<m_points.size();i++)
        for(int j=i+1;j<m_points.size();j++) {
            if (m_points[i].x >= m_points[j].x)
                std::swap(m_points[i], m_points[j]);
        }
    for(int i=0;i<m_points.size();i++)
        for(int j=i+1;j<m_points.size();j++) {
            if (m_points[i] == m_points[j]) {
                m_points.erase(m_points.begin() + j);
                i--;break;
            }
        }
}
void figure::setAccuracy(double accuracy){m_accuracy = accuracy;}
Point figure::getCenter(){return m_center;}
std::vector<Point> figure::getVertices(){return m_vertices;}
std::string figure::getType(){return m_type;}

triangle::triangle(Point A, Point B, Point C)
{
    m_type = "triangle";
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
void circle::setFrequency(double frequency){
    m_frequency = abs(frequency);
    m_points.clear();
    double tmp;
    for(double i = m_center.x - m_radius;i<=m_center.x+m_radius;i+=m_frequency)
    {
        tmp = sqrt(m_radius*m_radius - pow((i-m_center.x),2));
        m_points.push_back({i,m_center.y - tmp});
        m_points.push_back({i,m_center.y+tmp});
    }
    m_points.push_back({m_center.x+m_radius * cos(0),m_center.y+m_radius * sin(0)});
    m_points.push_back({m_center.x+m_radius * cos(acos(-1)/2),m_center.y+m_radius * sin(acos(-1)/2)});
    m_points.push_back({m_center.x+m_radius * cos(acos(-1)),m_center.y+m_radius * sin(acos(-1))});
    m_points.push_back({m_center.x+m_radius * cos(3*acos(-1)/2),m_center.y+m_radius * sin(3*acos(-1)/2)});
    for(int i=0;i<m_points.size();i++)
        for(int j=i+1;j<m_points.size();j++) {
            if (m_points[i].x == m_points[j].x && m_points[i].y == m_points[j].y) {
                m_points.erase(m_points.begin() + j);
                i--;break;
            }
        }
}
double circle::getSquare(){return acos(-1)*m_radius*m_radius;}
double circle::getPerimeter(){
    return 2*acos(-1)*m_radius;
}
double getSquare(Point A,Point B,Point C){
    double AB = sqrt(pow((A.x-B.x),2)+pow((A.y-B.y),2));
    double BC = sqrt(pow((B.x-C.x),2)+pow((B.y-C.y),2));
    double AC = sqrt(pow((A.x-C.x),2)+pow((A.y-C.y),2));
    double p = (AB+BC+AC)/2;
    return sqrt(p*(p-AB) * (p-BC) *(p-AC));}
double getSquare(Point A,Point B,Point C,Point D){
    segment tmpAB{A,B},tmpBC{B,C},tmpCD{C,D},tmpAD{A,D};
    double AB = tmpAB.getLength(), BC = tmpBC.getLength(),CD = tmpCD.getLength(),AD = tmpAD.getLength();
    double p = (AB+BC+CD+AD)/2;
    double cosAngle = cos((tmpAB.getAngle(&tmpBC)+tmpCD.getAngle(&tmpAD))/2);
    return sqrt((p-AB) * (p-BC) *(p-CD)*(p-AD)-AB * BC*CD*AD* pow(cosAngle,2));}

std::vector<Point> getIntersectionPoints(triangle* first, triangle* second)
{
    std::vector<Point> inPoints;
    if(first->getType()!="circle" && second->getType()!="circle") {
        for (auto &i: second->getBasePoint())
            for (auto &j: first->getBasePoint())
                if (abs(i.x - j.x) < basic_frequency * 10 && abs(i.y - j.y) < basic_frequency * 10)
                    inPoints.push_back(i);
    }
    else {
        for (auto &i: second->getBasePoint())
            for (auto &j: first->getBasePoint())
                if (abs(i.x - j.x) < basic_frequency && abs(i.y - j.y) < basic_frequency)
                    inPoints.push_back(i);
    }
    std::reverse(inPoints.begin(), inPoints.end());
    if(!inPoints.empty()) {
        int t=0;
        while(t<10) {
            for (int i = 0; i < inPoints.size() - 1; i++)
                for (int j = i + 1; j < inPoints.size(); j++) {
                    if (abs(inPoints[i].x - inPoints[j].x) < basic_frequency &&
                        abs(inPoints[i].y - inPoints[j].y) < basic_frequency) {
                        inPoints.erase(inPoints.begin() + j);
                        break;
                    }
                    if (inPoints.size()>2 && round(inPoints[i].x) == round(inPoints[j].x) &&
                        round(inPoints[i].y) == round(inPoints[j].y) ) {
                        int tmpX = round(inPoints[i].x);
                        int tmpY = round(inPoints[i].y);
                        if(abs(tmpX-inPoints[i].x)<=basic_frequency)
                            inPoints[i].x = round(inPoints[i].x);
                        if(abs(tmpY -  inPoints[i].y)<=basic_frequency)
                            inPoints[i].y = round(inPoints[i].y);
                        inPoints.erase(inPoints.begin() + j);
                        break;
                    }
                }
            t++;
        }
    }
    return  inPoints;
}
double getIntersectionSquare(triangle* first, triangle* second) {
    std::vector inPoints = getIntersectionPoints(first, second);
    if (inPoints.size() == 0)
        return 0;
    int tmp = inPoints.size();
    for (auto &i: second->getVertices()) {
        double tmp1 = getSquare(first->getVertices()[0],first->getVertices()[1],i);
        double tmp2 = getSquare(first->getVertices()[0],first->getVertices()[2],i);
        double tmp3 = getSquare(first->getVertices()[1],first->getVertices()[2],i);
        double a = first->getSquare();
        if(abs(tmp1+tmp2+tmp3-first->getSquare())<basic_frequency)
            inPoints.push_back(i);
    }
    if(inPoints.size()==tmp)
        for (auto &i: first->getVertices()) {
            double tmp1 = getSquare(second->getVertices()[0],second->getVertices()[1],i);
            double tmp2 = getSquare(second->getVertices()[0],second->getVertices()[2],i);
            double tmp3 = getSquare(second->getVertices()[1],second->getVertices()[2],i);
            double a = second->getSquare();
            double b = tmp1+ tmp2+tmp3;
            if(abs(tmp1+tmp2+tmp3-second->getSquare())<basic_frequency)
                inPoints.push_back(i);
        }
    for (int i = 0; i < inPoints.size() - 1; i++)
        for (int j = i + 1; j < inPoints.size(); j++) {
            if (abs(inPoints[i].x - inPoints[j].x) < basic_frequency &&
                abs(inPoints[i].y - inPoints[j].y) < basic_frequency) {
                inPoints.erase(inPoints.begin() + j);
                break;
            }}
    std::sort(inPoints.begin(), inPoints.end());

    if(inPoints.size()==3)
        return getSquare(inPoints[0],inPoints[1],inPoints[2]);
    else if(inPoints.size()==4)
        return getSquare(inPoints[0],inPoints[1],inPoints[2],inPoints[3]);
    else if(inPoints.size()==5)
        return getSquare(inPoints[0],inPoints[1],inPoints[2]) + getSquare(inPoints[1],inPoints[2],inPoints[3],inPoints[4]);
    else if(inPoints.size()==6)
        return getSquare(inPoints[0],inPoints[1],inPoints[2],inPoints[3]) + getSquare(inPoints[2],inPoints[3],inPoints[4],inPoints[5]);

}
