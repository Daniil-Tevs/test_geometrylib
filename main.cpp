#include "geometry.h"
#define  d_width -100
#define d_height 100
int main()
{
    quadrilateral test ({1,1},{2,3},{3,3},{4,1});
    for(int i=0;i<test.getGraph().size();i++)
       std::cout << test.getGraph()[i].x << " "<<test.getGraph()[i].y<<std::endl;
    test.setFrequency(0.5);
    for(int i=0;i<test.getGraph().size();i++)
        std::cout << test.getGraph()[i].x << " "<<test.getGraph()[i].y<<std::endl;
    std::cout<<test.getSquare();

}