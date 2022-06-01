#include "geometry.h"


bool Point::operator==(const Point& A) const {return (A.x == this->x) && (A.y == this->y);}
bool Point::operator>(const Point& A) const  {return (this->x>A.x);}
bool Point::operator<(const Point& A) const  {return (this->x < A.x);}

Point Point::operator+(const Point& A)
{
    this->x += A.x; this->y += A.y;
    return *this;
}
Point Point::operator-(const Point& A){
    this->x -= A.x; this->y -= A.y;
    return *this;
}

double getLength(Point A,Point B){return sqrt(pow((A.x-B.x),2)+pow((A.y-B.y),2));}

double getPerimeter(Point A,Point B,Point C)
{
    double AB = sqrt(pow((A.x-B.x),2)+pow((A.y-B.y),2));
    double BC = sqrt(pow((B.x-C.x),2)+pow((B.y-C.y),2));
    double AC = sqrt(pow((A.x-C.x),2)+pow((A.y-C.y),2));
    return AB+BC+AC;
}
double getPerimeter(Point A,Point B,Point C, Point D)
{
    double AB = getLength(A,B);
    double BC = getLength(B,C);
    double CD = getLength(C,D);
    double AD = getLength(A,D);
    return AB+BC+CD+AD;
}

double getSquare(Point A,Point B,Point C){
    double AB = sqrt(pow((A.x-B.x),2)+pow((A.y-B.y),2));
    double BC = sqrt(pow((B.x-C.x),2)+pow((B.y-C.y),2));
    double AC = sqrt(pow((A.x-C.x),2)+pow((A.y-C.y),2));
    double p = (AB+BC+AC)/2;
    return sqrt(p*(p-AB) * (p-BC) *(p-AC));}

double getSquare(Point A,Point B,Point C,Point D){
    std::vector<Point> tmp ={A,B,C,D};
    std::sort(tmp.begin(), tmp.end());
    segment tmpAB{tmp[0],tmp[1]},tmpBC{tmp[0],tmp[2]},tmpCD{tmp[1],tmp[3]},tmpAD{tmp[2],tmp[3]};
    double AB = tmpAB.getLength(), BC = tmpBC.getLength(),CD = tmpCD.getLength(),AD = tmpAD.getLength();
    double p = (AB+BC+CD+AD)/2;
    double cosAngle = cos((getAngle(&tmpAB,&tmpBC)+getAngle(&tmpCD,&tmpAD))/2);
    return sqrt((p-AB) * (p-BC) *(p-CD)*(p-AD)-AB * BC*CD*AD* pow(cosAngle,2));}

double getAngle(segment* AB,segment* CD)
{
    Point A = AB->getVertices().first; Point B = AB->getVertices().second;
    Point C = CD->getVertices().first; Point D = CD->getVertices().second;
    std::pair<double,double> firstV;
    std::pair<double,double> secondV;
    std::vector<Point> tmpV = getIntersectionPoints(AB,CD);
    double length = AB->getLength()*CD->getLength();
    if(tmpV.empty())
    {
        std::cout<<"don't intersection"<<std::endl;
        return 0;
    }
    else if(tmpV.size()>1)
        return 0;
    if(A.x == C.x && A.y == C.y)
    {
        firstV = std::make_pair(B.x-A.x,B.y-A.y);
        secondV = std::make_pair(D.x-C.x,D.y-C.y);
    }
    else if(B.x == D.x && B.y == D.y)
    {
        firstV = std::make_pair(A.x-B.x,A.y-B.y);
        secondV = std::make_pair(C.x-D.x,C.y-D.y);
    }
    else if(A.x == D.x && A.y == D.y)
    {
        firstV = std::make_pair(B.x-A.x,B.y-A.y);
        secondV = std::make_pair(C.x-D.x,C.y-D.y);
    }
    else if(B.x == C.x && B.y == C.y)
    {
        firstV = std::make_pair(A.x-B.x,A.y-B.y);
        secondV = std::make_pair(D.x-C.x,D.y-C.y);
    }
    else{
        if(tmpV.size()==1)
        {
            if(A.x == tmpV[0].x && A.y == tmpV[0].y) {length = getLength(B,tmpV[0]); firstV = std::make_pair(B.x - tmpV[0].x, B.y - tmpV[0].y); }
            else {length = getLength(A,tmpV[0]); firstV = std::make_pair(A.x - tmpV[0].x, A.y - tmpV[0].y); }
            if(C.x == tmpV[0].x && C.y == tmpV[0].y) {length*= getLength(D,tmpV[0]); secondV = std::make_pair(D.x - tmpV[0].x, D.y - tmpV[0].y); }
            else {length*= getLength(C,tmpV[0]); secondV = std::make_pair(C.x - tmpV[0].x, C.y - tmpV[0].y); }
        }
        else
            return 0;
    }
    double tmp = firstV.first*secondV.first + firstV.second*secondV.second;
    return acos(tmp/length);
}
double getAngleSegments(Point A, Point B, Point C, Point D)
{
    std::pair<double,double> firstV;
    std::pair<double,double> secondV;
    double length = getLength(A,B)*getLength(C,D);
    segment tmp_AB{A,B},tmp_CD{C,D};
    std::vector<Point> tmpV = getIntersectionPoints(&tmp_AB,&tmp_CD);
    if(tmpV.empty())
    {
        std::cout<<"don't intersection"<<std::endl;
        return 0;
    }
    else if(tmpV.size()>1)
        return 0;
    if(A.x == C.x && A.y == C.y)
    {
        firstV = std::make_pair(B.x-A.x,B.y-A.y);
        secondV = std::make_pair(D.x-C.x,D.y-C.y);
    }
    else if(B.x == D.x && B.y == D.y)
    {
        firstV = std::make_pair(A.x-B.x,A.y-B.y);
        secondV = std::make_pair(C.x-D.x,C.y-D.y);
    }
    else if(A.x == D.x && A.y == D.y)
    {
        firstV = std::make_pair(B.x-A.x,B.y-A.y);
        secondV = std::make_pair(C.x-D.x,C.y-D.y);
    }
    else if(B.x == C.x && B.y == C.y)
    {
        firstV = std::make_pair(A.x-B.x,A.y-B.y);
        secondV = std::make_pair(D.x-C.x,D.y-C.y);
    }
    else{
        if(tmpV.size()==1)
        {
            if(A.x == tmpV[0].x && A.y == tmpV[0].y) {length = getLength(B,tmpV[0]); firstV = std::make_pair(B.x - tmpV[0].x, B.y - tmpV[0].y); }
            else {length = getLength(A,tmpV[0]); firstV = std::make_pair(A.x - tmpV[0].x, A.y - tmpV[0].y); }
            if(C.x == tmpV[0].x && C.y == tmpV[0].y) {length*= getLength(D,tmpV[0]); secondV = std::make_pair(D.x - tmpV[0].x, D.y - tmpV[0].y); }
            else {length*= getLength(C,tmpV[0]); secondV = std::make_pair(C.x - tmpV[0].x, C.y - tmpV[0].y); }
        }
        else
            return 0;
    }
    double tmp = firstV.first*secondV.first + firstV.second*secondV.second;
    return acos(tmp/length);
}

std::vector<Point> getIntersectionPoints(segment* AB, segment* CD)
{
    std::vector<Point> inPoints;
    if(AB->getCollinear().first == 0 && CD->getCollinear().first!=0)
    {
        if(AB->getVertices().first.x<std::min(CD->getVertices().first.x,CD->getVertices().second.x) ||
           AB->getVertices().first.x>std::max(CD->getVertices().first.x,CD->getVertices().second.x))
            return inPoints;
        double tmpY = CD->functionY(AB->getVertices().first.x);
        double a = std::min(AB->getVertices().first.y,AB->getVertices().second.y);
        double b = std::max(AB->getVertices().first.y,AB->getVertices().second.y);
        if(tmpY>=a && tmpY<=b)
            inPoints.push_back({AB->getVertices().first.x,tmpY});
        return inPoints;
    }
    else if(CD->getCollinear().first == 0 && AB->getCollinear().first!=0)
    {
        if(CD->getVertices().first.x<std::min(AB->getVertices().first.x,AB->getVertices().second.x) ||
                CD->getVertices().first.x>std::max(AB->getVertices().first.x,AB->getVertices().second.x))
                return inPoints;
        double tmpY = AB->functionY(CD->getVertices().first.x);
        double a = std::min(CD->getVertices().first.y,CD->getVertices().second.y);
        double b = std::max(CD->getVertices().first.y,CD->getVertices().second.y);
        if(tmpY>=a && tmpY<=b)
            inPoints.push_back({CD->getVertices().first.x,tmpY});
        return inPoints;
    }
    else if(AB->getCollinear().first == 0 && CD->getCollinear().first==0)
    {
        if(AB->getVertices().first.x!=CD->getVertices().first.x)
            return inPoints;
        std::vector<Point> tmp;
        tmp.push_back(AB->getVertices().first); tmp.push_back(AB->getVertices().second);
        tmp.push_back(CD->getVertices().first); tmp.push_back(CD->getVertices().second);
        for(int i=0;i<tmp.size()-1;i++)
            for(int j = i+1;j<tmp.size();j++)
                if(tmp[i].y>tmp[j].y)
                    std::swap(tmp[i],tmp[j]);
        if(tmp[1]==tmp[2]) { inPoints.push_back(tmp[1]); return inPoints;}
        else if((tmp[0]==AB->getVertices().first && tmp[1]==AB->getVertices().second) || (tmp[0]==AB->getVertices().second && tmp[1]==AB->getVertices().first) ||
                (tmp[0]==CD->getVertices().first && tmp[1]==CD->getVertices().second) || (tmp[0]==CD->getVertices().second && tmp[1]==CD->getVertices().first))
        { return inPoints; }
        for(float y = tmp[1].y; y <= tmp[2].y; y+=basic_frequency) {
            inPoints.push_back({tmp[0].x, y}); }
        return inPoints;
    }
    std::vector<Point> tmp;
    tmp.push_back(AB->getVertices().first); tmp.push_back(AB->getVertices().second);
    tmp.push_back(CD->getVertices().first); tmp.push_back(CD->getVertices().second);
    std::sort(tmp.begin(), tmp.end());
    double a = tmp[1].x, b = tmp[2].x;
    for(float x = a; x <= b; x+=basic_frequency) {
        if (abs(AB->functionY(x) - CD->functionY(x)) < basic_accuracy) {
            float tmp_x = round(x),tmp_y = round(AB->functionY(x));
            if(abs(tmp_x-x)>basic_accuracy)
                tmp_x = x;
            if(abs(tmp_y-AB->functionY(x))>basic_accuracy)
                tmp_y = AB->functionY(x);
            inPoints.push_back({tmp_x, tmp_y}); }
    }
    if(inPoints.size()>=2 && inPoints.size()<5)
        while(inPoints.size()!=1)
            inPoints.erase(inPoints.begin());
    return inPoints;
}

std::vector<Point> getIntersectionPoints(segment* AB, circle* shape)
{
    std::vector<Point> inPoints;
    if(AB->getCollinear().first == 0)
    {
        double a = std::min(AB->getVertices().first.y,AB->getVertices().second.y);
        double b = std::max(AB->getVertices().first.y,AB->getVertices().second.y);
        std::pair<double,double>  tmpY = shape->function(AB->getVertices().first.x);
        if(tmpY.first>=a && tmpY.first<=b)
            inPoints.push_back({AB->getVertices().first.x,tmpY.first});
        else if(tmpY.second>=a && tmpY.second<=b)
            inPoints.push_back({AB->getVertices().first.x,tmpY.second});
        return inPoints;
    }
    std::vector<Point> tmp;
    double a = std::min(AB->getVertices().first.x,AB->getVertices().second.x);
    double b = std::max(AB->getVertices().first.x,AB->getVertices().second.x);
    for(float x = a; x <= b; x+=basic_frequency) {
        if ((abs(AB->functionY(x) - shape->function(x).first) < basic_accuracy) || (abs(AB->functionY(x) - shape->function(x).second) < basic_accuracy)) {
            float tmp_x = round(x),tmp_y = round(AB->functionY(x));
            if(abs(tmp_x-x)>basic_accuracy)
                tmp_x = x;
            if(abs(tmp_y-AB->functionY(x))>basic_accuracy)
                tmp_y = AB->functionY(x);
            inPoints.push_back({tmp_x, tmp_y}); }
    }
    if(inPoints.empty())
        return inPoints;
    for(int i=0;i<inPoints.size()-1;i++)
        for(int j=i+1;j<inPoints.size();j++)
            if(abs(inPoints[i].x-inPoints[j].x)<=basic_accuracy && abs(inPoints[i].y-inPoints[j].y)<=basic_accuracy) {
                inPoints.erase(inPoints.begin() + j);
                j--;
            }
    std::sort(inPoints.begin(), inPoints.end());
    return inPoints;
}

std::vector<Point> getIntersectionPoints(triangle* first, triangle* second)
{
    std::vector<Point> inPoints;
    for(auto& i : first->getSides())
        for(auto& j : second->getSides())
            for(auto& k : getIntersectionPoints(i,j))
                inPoints.push_back(k);
    if(inPoints.size()>12)
        return inPoints;
    for(int i=0;i<inPoints.size()-1;i++)
        for(int j=i+1;j<inPoints.size();j++)
            if(abs(inPoints[i].x-inPoints[j].x)<=basic_accuracy && abs(inPoints[i].y-inPoints[j].y)<=basic_accuracy) {
                inPoints.erase(inPoints.begin() + j);
                j--;
            }
    std::sort(inPoints.begin(), inPoints.end());
    return  inPoints;
}
std::vector<Point> getIntersectionPoints(quadrilateral* first, quadrilateral* second)
{
    std::vector<Point> inPoints;
    for(auto& i : first->getSides())
        for(auto& j : second->getSides())
            for(auto& k : getIntersectionPoints(i,j))
                inPoints.push_back(k);
    for(int i=0;i<inPoints.size()-1;i++)
        for(int j=i+1;j<inPoints.size();j++) {
            if (abs(inPoints[i].x - inPoints[j].x) <= basic_frequency*10 &&
                abs(inPoints[i].y - inPoints[j].y) <= basic_frequency*10 ) {
                inPoints.erase(inPoints.begin() + j);
                j--;
            }
        }
    std::sort(inPoints.begin(), inPoints.end());
    return  inPoints;
}
std::vector<Point> getIntersectionPoints(circle* first, circle* second)
{
    std::vector<Point> inPoints;
    std::vector<double> tmp;
    tmp.push_back(first->getRadius()*cos(acos(-1))+first->getCenter().x);
    tmp.push_back(first->getRadius()*cos(0)+first->getCenter().x);
    tmp.push_back(second->getRadius()*cos(acos(-1))+second->getCenter().x);
    tmp.push_back(second->getRadius()*cos(0)+second->getCenter().x);
    std::sort(tmp.begin(), tmp.end());
    double a=tmp[1], b = tmp[2];
    for(double x = a;x<=b;x+=basic_frequency)
    {
        std::pair<double,double> tmp1 = first->function(x);
        std::pair<double,double> tmp2 = second->function(x);
        if(abs(tmp1.first - tmp2.first)<basic_accuracy || abs(tmp1.second - tmp2.first)<basic_accuracy)
            inPoints.push_back({x, tmp2.first});
        if(abs(tmp1.first - tmp2.second)<basic_accuracy || abs(tmp1.second - tmp2.second)<basic_accuracy)
            inPoints.push_back({x,tmp2.second});
    }
    for(int i=0;i<inPoints.size()-1;i++)
        for(int j=i+1;j<inPoints.size();j++)
            if(abs(inPoints[i].x-inPoints[j].x)<=basic_accuracy && abs(inPoints[i].y-inPoints[j].y)<=basic_accuracy) {
                inPoints.erase(inPoints.begin() + j);
                j--;
            }
    std::sort(inPoints.begin(), inPoints.end());
    return  inPoints;
}
std::vector<Point> getIntersectionPoints(triangle* first, quadrilateral* second)
{
    std::vector<Point> inPoints;
    for(auto& i : first->getSides())
        for(auto& j : second->getSides())
            for(auto& k : getIntersectionPoints(i,j))
                inPoints.push_back(k);
    if(inPoints.empty())
        return inPoints;
    for(int i=0;i<inPoints.size()-1;i++)
        for(int j=i+1;j<inPoints.size();j++) {
            if (abs(inPoints[i].x - inPoints[j].x) <= basic_accuracy &&
                abs(inPoints[i].y - inPoints[j].y) <= basic_accuracy ) {
                inPoints.erase(inPoints.begin() + j);
                j--;
            }
        }
    std::sort(inPoints.begin(), inPoints.end());
    return  inPoints;
}
std::vector<Point> getIntersectionPoints(triangle* first, circle* second)
{
    std::vector<Point> inPoints;
    for(auto& i : first->getSides())
        for(auto& k : getIntersectionPoints(i,second))
            inPoints.push_back(k);
    if(inPoints.empty())
        return inPoints;
    for(int i=0;i<inPoints.size()-1;i++)
        for(int j=i+1;j<inPoints.size();j++)
            if(abs(inPoints[i].x-inPoints[j].x)<=basic_accuracy && abs(inPoints[i].y-inPoints[j].y)<=basic_accuracy) {
                inPoints.erase(inPoints.begin() + j);
                j--;
            }
    std::sort(inPoints.begin(), inPoints.end());
    return inPoints;
}

double getIntersectionSquare(triangle* first, triangle* second) {
    std::vector<Point> inPoints = getIntersectionPoints(first, second);
    int tmp = inPoints.size();

    for (auto &i: second->getVertices()) {
        if(std::find(inPoints.begin(), inPoints.end(), i) == inPoints.end()) {
            double tmp1 = getSquare(first->getVertices()[0], first->getVertices()[1], i);
            double tmp2 = getSquare(first->getVertices()[0], first->getVertices()[2], i);
            double tmp3 = getSquare(first->getVertices()[1], first->getVertices()[2], i);
            if (abs(tmp1 + tmp2 + tmp3 - first->getSquare()) < basic_frequency)
                inPoints.push_back(i);
        }
    }
    if(inPoints.size()==tmp)
        for (auto &i: first->getVertices()) {
            if(std::find(inPoints.begin(), inPoints.end(), i) == inPoints.end()) {
                double tmp1 = getSquare(second->getVertices()[0], second->getVertices()[1], i);
                double tmp2 = getSquare(second->getVertices()[0], second->getVertices()[2], i);
                double tmp3 = getSquare(second->getVertices()[1], second->getVertices()[2], i);
                if (abs(tmp1 + tmp2 + tmp3 - second->getSquare()) < basic_frequency)
                    inPoints.push_back(i);
            }
        }

    std::sort(inPoints.begin(), inPoints.end());
    for(int i=0;i<inPoints.size()-1;i++)
        for(int j=i+1;j<inPoints.size();j++)
        {
            if(inPoints[i].x == inPoints[j].x)
                if(inPoints[i].y>inPoints[j].y)
                    std::swap(inPoints[i],inPoints[j]);
        }
    if(inPoints.size()==3)
        return getSquare(inPoints[0],inPoints[1],inPoints[2]);
    else if(inPoints.size()==4)
        return getSquare(inPoints[0],inPoints[1],inPoints[3],inPoints[2]);
    else if(inPoints.size()==5)
        return getSquare(inPoints[0],inPoints[1],inPoints[2]) + getSquare(inPoints[1],inPoints[2],inPoints[4],inPoints[3]);
    else if(inPoints.size()==6)
        return getSquare(inPoints[0],inPoints[1],inPoints[3],inPoints[2]) + getSquare(inPoints[2],inPoints[3],inPoints[5],inPoints[4]);
    return 0;
}
double getIntersectionSquare(quadrilateral* first, quadrilateral* second) {
    std::vector<Point> inPoints = getIntersectionPoints(first, second);
    int tmp = inPoints.size();

    for (auto &i: second->getVertices()) {
        if(std::find(inPoints.begin(), inPoints.end(), i) == inPoints.end()) {
            double tmp1 = getSquare(first->getVertices()[0], first->getVertices()[1], i);
            double tmp2 = getSquare(first->getVertices()[1], first->getVertices()[2], i);
            double tmp3 = getSquare(first->getVertices()[2], first->getVertices()[3], i);
            double tmp4 = getSquare(first->getVertices()[3], first->getVertices()[0], i);
            if (abs(tmp1 + tmp2 + tmp3 + tmp4 - first->getSquare()) < basic_frequency)
                inPoints.push_back(i);
        }
    }
    for (auto &i: first->getVertices()) {
        if(std::find(inPoints.begin(), inPoints.end(), i) == inPoints.end()) {
            double tmp1 = getSquare(second->getVertices()[0], second->getVertices()[1], i);
            double tmp2 = getSquare(second->getVertices()[1], second->getVertices()[2], i);
            double tmp3 = getSquare(second->getVertices()[2], second->getVertices()[3], i);
            double tmp4 = getSquare(second->getVertices()[3], second->getVertices()[0], i);
            if (abs(tmp1 + tmp2 + tmp3 + tmp4 - second->getSquare()) < basic_frequency)
                inPoints.push_back(i);
        }
    }

    std::sort(inPoints.begin(), inPoints.end());
    for(int i=0;i<inPoints.size()-1;i++)
        for(int j=i+1;j<inPoints.size();j++)
        {
            if(inPoints[i].x == inPoints[j].x)
                if(inPoints[i].y>inPoints[j].y)
                    std::swap(inPoints[i],inPoints[j]);
        }

    if(inPoints.size()==3)
        return getSquare(inPoints[0],inPoints[1],inPoints[2]);
    else if(inPoints.size()==4)
        return getSquare(inPoints[0],inPoints[1],inPoints[2],inPoints[3]);
    else if(inPoints.size()==5)
        return getSquare(inPoints[0],inPoints[1],inPoints[2]) + getSquare(inPoints[1],inPoints[2],inPoints[4],inPoints[3]);
    else if(inPoints.size()==6)
        return getSquare(inPoints[0],inPoints[1],inPoints[2],inPoints[3]) + getSquare(inPoints[2],inPoints[3],inPoints[4],inPoints[5]);
    else if(inPoints.size() == 7)
    {}
    else if(inPoints.size() == 8)
    {

    }
    return 0;
}
double getIntersectionSquare(circle* first, circle* second) {
    std::vector<Point> inPoints = getIntersectionPoints(first, second);
    if(inPoints.size() <= 1)
    {
        if(first->getCenter().x > second->getRadius()*cos(acos(-1))+second->getCenter().x &&
            first->getCenter().x < second->getRadius()+second->getCenter().x)
            return first->getSquare();
        else if(second->getCenter().x > first->getRadius()*cos(acos(-1))+first->getCenter().x &&
                second->getCenter().x < first->getRadius()+first->getCenter().x)
            return second->getSquare();
        else
            return 0;
    }
    else if(inPoints.size()==2) {
        double centers_and_inPoints_square = getSquare(first->getCenter(),inPoints[0], second->getCenter(), inPoints[1]);
        segment sec1_AB{first->getCenter(), inPoints[0]};
        segment sec1_CD{first->getCenter(), inPoints[1]};
        double square_sec1 = getAngle(&sec1_AB, &sec1_CD) / 2 * pow(first->getRadius(), 2);
        segment sec2_AB{second->getCenter(), inPoints[0]};
        segment sec2_CD{second->getCenter(), inPoints[1]};
        double square_sec2 = getAngle(&sec2_AB, &sec2_CD) / 2 * pow(second->getRadius(), 2);
        if(square_sec1>first->getSquare()/2)
            square_sec1 = first->getSquare()-square_sec1;
        if(square_sec2>second->getSquare()/2)
            square_sec2 = second->getSquare()-square_sec2;
        if(abs(second->function(first->getCenter().x).first-first->getCenter().y)<basic_accuracy || abs(second->function(first->getCenter().x).second-first->getCenter().y)<basic_accuracy)
        {
            double tmp = square_sec2-centers_and_inPoints_square;
            return square_sec1+tmp;
        }
        else if(abs(first->function(second->getCenter().x).first-second->getCenter().y)<basic_accuracy || abs(first->function(second->getCenter().x).second-second->getCenter().y)<basic_accuracy)
        {
            double tmp = square_sec1-centers_and_inPoints_square;
            return square_sec2+tmp;
        }
        else if ((first->getCenter().x > second->getRadius() * cos(acos(-1)) + second->getCenter().x &&
             first->getCenter().x < second->getRadius() + second->getCenter().x)){
            double tmp = square_sec2 - centers_and_inPoints_square;
            return first->getSquare()-(square_sec1 - tmp);
        }
        else if ((second->getCenter().x > first->getRadius() * cos(acos(-1)) + first->getCenter().x &&
                     second->getCenter().x < first->getRadius() + first->getCenter().x)) {
            double tmp = square_sec1 - centers_and_inPoints_square;
            return second->getSquare()-(square_sec2 - tmp);
        }
        else
        {
            double tmp = centers_and_inPoints_square - square_sec1;
            return square_sec2 - tmp;
        }
    }
    else
        return first->getSquare();

}
double getIntersectionSquare(triangle* first, quadrilateral* second) {
    auto * tmp1 = new triangle{second->getVertices()[0],second->getVertices()[1],second->getVertices()[2]};
    auto * tmp2 = new triangle{second->getVertices()[1],second->getVertices()[2],second->getVertices()[3]};
    double square1 = getIntersectionSquare(first,tmp1);
    double square2 = getIntersectionSquare(first,tmp2);
    delete tmp1;
    delete tmp2;
    return square1+square2;
}

double getIntersectionSquare(triangle* first, circle* second) {
    std::vector<Point> inPoints = getIntersectionPoints(first, second);
    double X, Y;
    if(inPoints.size()<=1 || inPoints.size()==3)
    {
        X = first->getVertices()[0].x;
        Y = first->getVertices()[0].y;
        double a = second->getCenter().x-second->getRadius();
        double b = second->getCenter().x+second->getRadius();
        if(X>=second->getCenter().x-second->getRadius() && X<=second->getCenter().x+second->getRadius())
            if(Y>=second->getCenter().y-second->getRadius() && Y<=second->getCenter().y+second->getRadius())
                return first->getSquare();
        Point C{second->getCenter().x,second->getCenter().y};
        double tmp1 = getSquare(C, first->getVertices()[0], first->getVertices()[1]);
        double tmp2 = getSquare(C, first->getVertices()[0], first->getVertices()[2]);
        double tmp3 = getSquare(C, first->getVertices()[1], first->getVertices()[2]);
        if (abs(tmp1 + tmp2 + tmp3 - first->getSquare()) < basic_frequency)
            return second->getSquare();
        return 0;
    }
    else if(inPoints.size() == 2)
    {
        segment AB{inPoints[0],second->getCenter()};
        segment CD{inPoints[1],second->getCenter()};
        double tmpSquare = pow(second->getRadius(),2)* getAngle(&AB,&CD)/2 - getSquare(inPoints[0],inPoints[1],second->getCenter());
        for(auto& i : first->getVertices())
        {
            if(std::find(inPoints.begin(), inPoints.end(), i)==inPoints.end()) {
                X = i.x;
                Y = i.y;
                if (X >= second->getRadius() - second->getCenter().x &&
                    X <= second->getRadius() + second->getCenter().x)
                    if (Y >= second->getRadius() - second->getCenter().y &&
                        Y <= second->getRadius() + second->getCenter().y)
                        inPoints.push_back(i);
            }
        }
        if(inPoints.size() == 3 && std::find(first->getVertices().begin(), first->getVertices().end(), inPoints[0])!=first->getVertices().end()
        && std::find(first->getVertices().begin(), first->getVertices().end(), inPoints[1])!=first->getVertices().end())
            return first->getSquare();
        else if(inPoints.size() == 2 && std::find(first->getVertices().begin(), first->getVertices().end(), inPoints[0])!=first->getVertices().end()
           && std::find(first->getVertices().begin(), first->getVertices().end(), inPoints[1])!=first->getVertices().end())
            return tmpSquare;
        if(inPoints.size()==3)
            return getSquare(inPoints[0],inPoints[1],inPoints[2])+tmpSquare;
        else if(inPoints.size()==4)
            return getSquare(inPoints[0],inPoints[1],inPoints[2],inPoints[3])+tmpSquare;
        return 0;
    }
    std::sort(inPoints.begin(), inPoints.end());
    for(int i=0;i<inPoints.size()-1;i++)
        for(int j=i+1;j<inPoints.size();j++)
        {
            if(inPoints[i].x == inPoints[j].x)
                if(inPoints[i].y>inPoints[j].y)
                    std::swap(inPoints[i],inPoints[j]);
        }
    if(inPoints.size()==4)
        return getSquare(inPoints[0],inPoints[1],inPoints[3],inPoints[2]);
    else if(inPoints.size()==5)
        return getSquare(inPoints[0],inPoints[1],inPoints[2]) + getSquare(inPoints[1],inPoints[2],inPoints[4],inPoints[3]);
    else if(inPoints.size()==6)
        return getSquare(inPoints[0],inPoints[1],inPoints[3],inPoints[2]) + getSquare(inPoints[2],inPoints[3],inPoints[5],inPoints[4]);
    return 0;
}

double getIntersectionSquare(quadrilateral* first, circle* second) {
    auto * tmp1 = new triangle{first->getVertices()[0],first->getVertices()[1],first->getVertices()[2]};
    auto * tmp2 = new triangle{first->getVertices()[1],first->getVertices()[2],first->getVertices()[3]};
    double square1 = getIntersectionSquare(tmp1,second);
    double square2 = getIntersectionSquare(tmp2,second);
    delete tmp1;
    delete tmp2;
    return square1+square2;
}