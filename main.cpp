#include "geometry.h"

int main()
{
    triangle test ({0,4},{-3,-2},{3,-2});
    triangle r({-4,2},{4,2},{0,-2});
    segment AB{{0,6},{0,2}}, CD{{0,0},{0,3}};

    std::cout<<getIntersectionSquare(&test,&r);
}