#include "geometry.h"

int main()
{
//    triangle test ({-2,0},{0,3},{2,0});
//    triangle r ({-3,2},{0,-1},{3,2});
////    triangle r({-4,2},{4,2},{0,-2});
//    segment AB{{0,1},{-1,-2}}, CD{{0,0},{10,0}};
//    std::cout<<getAngle(&AB,&CD)*180/acos(-1);

//    quadrilateral test ({0,0},{0,2},{2,2},{2,0});
//    quadrilateral r({1,1},{1,3},{3,3},{3,1});
//    for(auto& i : getIntersectionPoints(&test,&r))
//        std::cout<<i<<std::endl;


//    circle test ({0,0},1);
////    circle r({0,1.5},2);
//    for(auto& i : getIntersectionPoints(&test,&r))
//        std::cout<<i<<std::endl;
//    std::cout<<getIntersectionSquare(&test,&r);
//    std::cout<<getIntersectionSquare(&test,&r);
//    std::cout<<test.getPerimeter()<<std::endl;
//    std::cout<<test.getSquare()<<std::endl;
//    std::cout<<test.getCenter()<<std::endl;
    quadrilateral test{{-2,0},{0,2},{2,0},{0,-2}};
//    quadrilateral r{{-1,-1},{-1,2},{3,2},{3,-2}};
    circle r({0,0},1);
    std::cout<<getIntersectionSquare(&test,&r);
}