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
    segment tmpAB{A,B},tmpBC{B,C},tmpCD{C,D},tmpAD{A,D};
    double AB = tmpAB.getLength(), BC = tmpBC.getLength(),CD = tmpCD.getLength(),AD = tmpAD.getLength();
    double p = (AB+BC+CD+AD)/2;
    double cosAngle = cos((tmpAB.getAngle(&tmpBC)+tmpCD.getAngle(&tmpAD))/2);
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
        double tmpY = CD->functionY(AB->getVertices().first.x);
        double a = std::min(AB->getVertices().first.y,AB->getVertices().second.y);
        double b = std::max(AB->getVertices().first.y,AB->getVertices().second.y);
        if(tmpY>=a && tmpY<=b)
            inPoints.push_back({AB->getVertices().first.x,tmpY});
        return inPoints;
    }
    else if(CD->getCollinear().first == 0 && AB->getCollinear().first!=0)
    {
        double tmpY = AB->functionY(CD->getVertices().first.x);
        double a = std::min(CD->getVertices().first.y,CD->getVertices().second.y);
        double b = std::max(CD->getVertices().first.y,CD->getVertices().second.y);
        if(tmpY>=a && tmpY<=b)
            inPoints.push_back({CD->getVertices().first.x,tmpY});
        return inPoints;
    }
    else if(AB->getCollinear().first == 0 && CD->getCollinear().first==0)
    {
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

//std::vector<Point> getIntersectionPoints(figure* first, figure* second)
//{
//    std::vector<Point> inPoints;
//    if(first->getType()!="circle" && second->getType()!="circle") {
//        for (auto &i: second->getPoints())
//            for (auto &j: first->getPoints())
//                if (abs(i.x - j.x) < basic_frequency * 10 && abs(i.y - j.y) < basic_frequency * 10)
//                    inPoints.push_back(i);
//    }
//    else {
//        for (auto &i: second->getPoints())
//            for (auto &j: first->getPoints())
//                if (abs(i.x - j.x) < basic_frequency && abs(i.y - j.y) < basic_frequency)
//                    inPoints.push_back(i);
//    }
//    std::reverse(inPoints.begin(), inPoints.end());
//    if(!inPoints.empty()) {
//        int t=0;
//        while(t<10) {
//            for (int i = 0; i < inPoints.size() - 1; i++)
//                for (int j = i + 1; j < inPoints.size(); j++) {
//                    if (abs(inPoints[i].x - inPoints[j].x) < basic_frequency &&
//                        abs(inPoints[i].y - inPoints[j].y) < basic_frequency) {
//                        inPoints.erase(inPoints.begin() + j);
//                        break;
//                    }
//                    if (inPoints.size()>2 && round(inPoints[i].x) == round(inPoints[j].x) &&
//                        round(inPoints[i].y) == round(inPoints[j].y) ) {
//                        int tmpX = round(inPoints[i].x);
//                        int tmpY = round(inPoints[i].y);
//                        if(abs(tmpX-inPoints[i].x)<=basic_frequency)
//                            inPoints[i].x = round(inPoints[i].x);
//                        if(abs(tmpY -  inPoints[i].y)<=basic_frequency)
//                            inPoints[i].y = round(inPoints[i].y);
//                        inPoints.erase(inPoints.begin() + j);
//                        break;
//                    }
//                }
//            t++;
//        }
//    }
//    return  inPoints;
//}
//
double getIntersectionSquare(triangle* first, triangle* second) {
    std::vector inPoints = getIntersectionPoints(first, second);
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
        return getSquare(inPoints[0],inPoints[1],inPoints[2],inPoints[3]);
    else if(inPoints.size()==5)
        return getSquare(inPoints[0],inPoints[1],inPoints[2]) + getSquare(inPoints[1],inPoints[2],inPoints[4],inPoints[3]);
    else if(inPoints.size()==6)
        return getSquare(inPoints[0],inPoints[1],inPoints[2],inPoints[3]) + getSquare(inPoints[2],inPoints[3],inPoints[4],inPoints[5]);

}
//std::vector<Point> getIntersectionPoints(segment* first, segment* segment)
//{
//    std::vector<Point> inPoints;
//    for (auto &i: first->getGraph())
//        for (auto &j: segment->getGraph())
//            if (abs(i.x - j.x) < basic_accuracy && abs(i.y - j.y) < basic_accuracy)
//                inPoints.push_back(i);
//    std::reverse(inPoints.begin(), inPoints.end());
//    if(!inPoints.empty()) {
//        int t=0;
//        while(t<2) {
//            for (int i = 0; i < inPoints.size() - 1; i++)
//                for (int j = i + 1; j < inPoints.size(); j++) {
//                    if (abs(inPoints[i].x - inPoints[j].x) < basic_accuracy &&
//                        abs(inPoints[i].y - inPoints[j].y) < basic_accuracy) {
//                        inPoints.erase(inPoints.begin() + j);
//                        break;
//                    }
//                    if (round(inPoints[i].x) == round(inPoints[j].x) &&
//                        round(inPoints[i].y) == round(inPoints[j].y)) {
//                        int tmpX = round(inPoints[i].x);
//                        int tmpY = round(inPoints[i].y);
//                        if(abs(tmpX-inPoints[i].x)<=basic_accuracy*5)
//                            inPoints[i].x = round(inPoints[i].x);
//                        if(abs(tmpY -  inPoints[i].y)<=basic_accuracy*5)
//                            inPoints[i].y = round(inPoints[i].y);
//                        inPoints.erase(inPoints.begin() + j);
//                        break;
//                    }
//                }
//            t++;
//        }
//    }
//    return inPoints;
//}

//    if(!inPoints.empty()) {
//        int t=0;
//        while(t<10) {
//            for (int i = 0; i < inPoints.size() - 1; i++)
//                for (int j = i + 1; j < inPoints.size(); j++) {
//                    if (abs(inPoints[i].x - inPoints[j].x) < basic_frequency &&
//                        abs(inPoints[i].y - inPoints[j].y) < basic_frequency) {
//                        inPoints.erase(inPoints.begin() + j);
//                        break;
//                    }
//                    if (inPoints.size()>2 && round(inPoints[i].x) == round(inPoints[j].x) &&
//                        round(inPoints[i].y) == round(inPoints[j].y) ) {
//                        int tmpX = round(inPoints[i].x);
//                        int tmpY = round(inPoints[i].y);
//                        if(abs(tmpX-inPoints[i].x)<=basic_frequency)
//                            inPoints[i].x = round(inPoints[i].x);
//                        if(abs(tmpY -  inPoints[i].y)<=basic_frequency)
//                            inPoints[i].y = round(inPoints[i].y);
//                        inPoints.erase(inPoints.begin() + j);
//                        break;
//                    }
//                }
//            t++;
//        }
//    }