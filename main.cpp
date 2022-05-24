#include "geometry.h"
int main()
{
    quadrilateral test ({-2,-2},{-2,2},{2,2},{2,-2});
    for(int i=0;i<test.getGraph().size();i++)
       std::cout << test.getGraph()[i].x << " "<<test.getGraph()[i].y<<std::endl;
//    test.setFrequency(0.5);
//    for(int i=0;i<test.getGraph().size();i++)
//        std::cout << test.getGraph()[i].x << " "<<test.getGraph()[i].y<<std::endl;
    std::cout<<test.getPerimeter();

    circle test2({0,0},1);
    std::cout<<test2.getPerimeter();

}