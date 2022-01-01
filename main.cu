#include <vector>
#include <iostream>

__global__ class AABB{
public:
    double3 min;
    double3 max;

    double3 max1;
    double3 max2;
    double3 max3;
    double3 max4;
    double3 max5;

    int3 asd;
};

void bla1(std::vector<AABB> &a)
{
    for(int i = 0; i < 10; ++i)
        a.emplace_back();
}

void bla(std::vector<AABB> &x)
{
    bla1(x);
    bla1(x);
    bla1(x);
}





int main()
{
    std::vector<AABB> boxes;
    bla(boxes);
    std::cout<<boxes.size()<<std::endl;
    return 0;
}