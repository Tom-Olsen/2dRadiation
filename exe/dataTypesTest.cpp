#include <iostream>
#include "../src/DataTypes.hh"
using namespace std;



void TestCoord()
{
    Coord xy(1,1);
    xy.Print("xy");
    cout << "xy:" << xy << endl;
    PrintDouble(xy[1],"x");
    PrintDouble(xy[2],"y");
    cout << endl;
}
void TestTensor2()
{
    cout << "Tensor2" << endl;
    Tensor2 a(1,2);
    a.Print("a");
    cout << "a:" << a << endl;;
    for(int i=1; i<3; i++)
        PrintDouble(a[i],"a[" + std::to_string(i) + "]");
    Tensor2 b(4); 
    b.Print("b");
    cout << "b:" << b << endl;;
    for(int i=1; i<3; i++)
        PrintDouble(b[i],"b[" + std::to_string(i) + "]");
    cout << endl;
}
void TestTensor3()
{
    cout << "Tensor3" << endl;
    Tensor3 a(1,2,3);
    cout << "a:" << a << endl;;
    a.Print("a4");
    for(int i=0; i<3; i++)
        PrintDouble(a[i],"a[" + std::to_string(i) + "]");
    Tensor3 b(4); 
    cout << "b:" << b << endl;;
    b.Print("b");
    for(int i=0; i<3; i++)
        PrintDouble(b[i],"b[" + std::to_string(i) + "]");
    cout << endl;
}
void TestTensor2x2()
{
    cout << "Tensor2x2" << endl;
    Tensor2x2 A(0,1, 2,3);
    A.Print("A");
    cout << "A:" << A << endl;;
    for(int i=1; i<3; i++)
        for(int j=1; j<3; j++)
            PrintDouble(A[{i,j}],"A[" + std::to_string(i) + "," + std::to_string(j) + "]");
    Tensor2x2 B(9);
    B.Print("B");
    cout << "B:" << B << endl;;
    for(int i=1; i<3; i++)
        for(int j=1; j<3; j++)
            PrintDouble(B[{i,j}],"B[" + std::to_string(i) + "," + std::to_string(j) + "]");
    cout << endl;
}
void TestTensor3x3()
{
    cout << "Tensor3x3" << endl;
    Tensor3x3 A(0,1,2, 3,4,5, 6,7,8);
    A.Print("A");
    cout << "A:" << A << endl;;
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            PrintDouble(A[{i,j}],"A[" + std::to_string(i) + "," + std::to_string(j) + "]");
    Tensor3x3 B(1.234);
    B.Print("B");
    cout << "B:" << B << endl;;
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            PrintDouble(B[{i,j}],"B[" + std::to_string(i) + "," + std::to_string(j) + "]");
    cout << endl;
}
void TestTensor2x2x2()
{
    cout << "Tensor2x2x2" << endl;
    Tensor2x2x2 A(0);
    for(int i=1; i<3; i++)
        for(int j=1; j<3; j++)
            for(int k=1; k<3; k++)
                A[{i,j,k}] = (k-1) + 2*(j-1) + 4*(i-1);
    A.Print("A");
    cout << "A:" << A << endl;;
    for(int i=1; i<3; i++)
        for(int j=1; j<3; j++)
            for(int k=1; k<3; k++)
                PrintDouble(A[{i,j,k}],"A[" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + "]");
    cout << endl;
}
void TestTensor3x3x3()
{
    cout << "Tensor3x3x3" << endl;
    Tensor3x3x3 A(0);
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            for(int k=0; k<3; k++)
                A[{i,j,k}] = k + 3*j + 9*i;
    A.Print("A");
    cout << "A:" << A << endl;;
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            for(int k=0; k<3; k++)
                PrintDouble(A[{i,j,k}],"A[" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + "]");
    cout << endl;
}
void TestRotationMatrix()
{
    RotationMatrix rot(M_PI/4.0);
    Tensor2 v(1.0,0.0);

    v.Print("v");
    (rot * v).Print("vRot");
}



int main()
{
    TestCoord();
    TestTensor2();
    TestTensor3();
    TestTensor2x2();
    TestTensor3x3();
    TestTensor2x2x2();
    TestTensor3x3x3();
    TestRotationMatrix();
}