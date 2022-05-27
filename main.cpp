#include "geometry.h"

int main()
{
    triangle test ({0,0},{1,3},{-1,3});
    triangle r({-1,1},{0,4},{1,1});
//    circle test2 ({0,0},2);
//    circle test ({2,0},2);
////    for(int i=0;i<test.getGraph().size();i++)
////        std::cout << test.getGraph()[i].x << " "<<test.getGraph()[i].y<<std::endl;
////    std::cout<<"tr 2 ***********"<<std::endl;
////    for(int i=0;i<test2.getGraph().size();i++)
////       std::cout << test2.getGraph()[i].x << " "<<test2.getGraph()[i].y<<std::endl;
////    test.setFrequency(0.5);
////    for(int i=0;i<test.getGraph().size();i++)
////        std::cout << test.getGraph()[i].x << " "<<test.getGraph()[i].y<<std::endl;
////    std::cout<<test.getPerimeter();
//    std::vector<Point> tmp;
//    tmp = getIntersectionPoints(&test,&r);
//    for(auto& i : tmp)
//        std::cout<<i<<std::endl;
//    std::cout<<round(-2.121);
    segment AB{{-1,-1},{1,1}}, CD{{-1,1},{0,0}};
//    std::cout<<getAngleSegments({-1,-1}, {1,1}, {-1,1},{0,0});
    std::cout<<getAngle(&AB,&CD);
}