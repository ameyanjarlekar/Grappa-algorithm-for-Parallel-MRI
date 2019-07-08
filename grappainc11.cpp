//load('C:\Users\lenovo\Desktop\Raw_data_recon.mat')
//ksp got by fft of image data
#include <iostream>
#include <vector>
#include <complex>
#include "imag_coil1.h"
#include "imag_coil2.h"
#include "real_coil1.h"
#include "real_coil2.h"

using namespace std;
/*class complex{
public:
    float x;
    float y;
    complex() {
        x = 0;
        y = 0;
    }
    complex (float a,float b)
    { x =a;
        y = b;}
    complex( complex &a)
    {x = a.x;
        y = a.y;}
    complex conjugate()
    {
    x = x;
    y = -y;
     return *this;}
     void print()
     {
       cout<<x<<" "<<"+"<<" "<<"i"<<y;
     return;}

    complex operator = (complex & a)
    {complex b(a.x,a.y);
        return b;}

};

    complex conjugate(complex &a)
    {complex b;
    b.x = a.x;
    b.y = -a.y;
     return b;}

ostream & operator << (ostream &ost,complex &a)
{ost <<a.x<<" "<<"+"<<" "<<"i"<<a.y<<" ";
return ost;}

istream & operator >> (istream &cin,complex &a)
{char b;
    cin >>a.x;
    cin>>b;
    cin>>b;
    cin>>a.y;
return cin;}


    complex operator + (complex & c,complex & a){
    complex b(c.x+a.x,c.x+a.y);
    //b.x = c.x+a.x;
    //b.y = c.x+a.y;
    return b;}
    complex operator - (complex & c,complex & a){
        complex b(c.x-a.x,c.x-a.y);
    //b.x = c.x-a.x;
    //b.y = c.x-a.y;
    return b;
    }    
        complex operator * (complex  c,complex  a){
        complex b(a.x * c.x- c.y* a.y,c.x*a.y + c.y*a.x);
    //b.x = a.x * c.x- c.y* a.y;
    //b.y = c.x*a.y + c.y*a.x;
    return b;
    }

        complex operator * (float  c,complex  a){
        complex b(a.x*c,a.y*c);
    //b.x = a.x * c.x- c.y* a.y;
    //b.y = c.x*a.y + c.y*a.x;
    return b;
    }

        complex operator / (complex & c,complex & a){
        complex b((c.x*a.x + c.y*b.y)/(a.x*a.x + a.y*a.y),(a.y*b.x - a.x*b.y)/(a.x*a.x + a.y*a.y));
        //b.x = (c.x*a.x + c.y*b.y)/(a.x*a.x + a.y*a.y);
        //b.y = (a.y*b.x - a.x*b.y)/(a.x*a.x + a.y*a.y);
    return b;
    }
*/
// calculate the cofactor of element (row,col)
 

void GJinverse (complex<float>  b[][12])
{
    int n=12;
    complex<float> a[12][24];
    //float a[10][10] = { 0 }, d;
/*    cout << "No of equations ? ";
    cin >> n;
    cout << "Read all coefficients of matrix with b matrix too " << endl;
    for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++)
            cin >> a[i][j];*/
 
    for (int i = 0; i < n; i++)
        for (int j = 0; j < 2 * n; j++)
            if (j == (i + n))
                a[i][j] = 1.0;
            else 
                a[i][j] = 0.0;

     for (int i = 0; i < n; i++)
        for (int j = 0; j <  n; j++)
            a[i][j] = b[i][j];

        for (int i = 0; i < n; i++)
            {for (int j = 0; j < 2 * n; j++) 
                {cout<<a[i][j]<<" ";}
            cout<<endl;}   
    /************** partial pivoting **************/
/*    for (int i = n-1; i > 0; i--)
    {
        if (a[i - 1][0] < a[i][0])
            for (int j = 0; j < n * 2; j++)
            {float d;
                d = a[i][j];
                a[i][j] = a[i - 1][j];
                a[i - 1][j] = d;
            }
    }
    cout << "pivoted output: " << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n * 2; j++)
            cout << a[i][j] << "    ";
        cout << endl;
    }*/
    /********** reducing to diagonal  matrix ***********/
 
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            if (j != i)
            {complex<float> d;
                d = a[j][i] / a[i][i];
                for (int k = 0; k < n * 2; k++)
                    a[j][k] -= a[i][k] * d;
            }
    }
    cout<<"diagonal matrix"<<endl;
        for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n * 2; j++)
            cout << a[i][j] << "    ";
        cout << endl;
    }
    /************** reducing to unit matrix *************/
    for (int i = 0; i < n; i++)
    {complex<float> d;
        d = a[i][i];
        for (int j = 0; j < n * 2; j++)
            a[i][j] = a[i][j] / d;
    }
 
    cout << "your solutions: " << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = n ; j < n * 2; j++)
            {b[i][j-n] = a[i][j];
              cout << a[i][j] << "    ";}
          
        cout << endl;
    }
return;
}


void getCofactor(complex<float> A[][12], complex<float> temp[][12], int p, int q, int n) 
{ 
    int i = 0, j = 0; 
  
    // Looping for each element of the matrix 
    for (int row = 0; row < n; row++) 
    { 
        for (int col = 0; col < n; col++) 
        { 
            //  Copying into temporary matrix only those element 
            //  which are not in given row and column 
            if (row != p && col != q) 
            { 
                temp[i][j++] = A[row][col]; 
  
                // Row is filled, so increase row index and 
                // reset col index 
                if (j == n - 1) 
                { 
                    j = 0; 
                    i++; 
                } 
            } 
        } 
    } 
} 

  complex<float> determinant(complex<float> A[][12], int n) 
{ 
    complex<float> D = 0; // Initialize result 
  
    //  Base case : if matrix contains single element 
    if (n == 1) 
        return A[0][0]; 
  
    complex<float> temp[12][12]; // To store cofactors 
  
    float sign = 1.0;  // To store sign multiplier 
  
     // Iterate for each element of first row 
    for (int f = 0; f < n; f++) 
    { 
        // Getting Cofactor of A[0][f] 
        getCofactor(A, temp, 0, f, n); 
        D = D + sign * A[0][f] * determinant(temp, n - 1); 
  
        // terms are to be added with alternate sign 
        sign = -sign; 
    } 
  
    return D; 
} 


void adjoint(complex<float> A[][12],complex<float> adj[][12]) 
{ 
/*    if (N == 1) 
    { 
        adj[0][0] = 1; 
        return; 
    } 
  */
    // temp is used to store cofactors of A[][] 
    float sign = 1.0;
    complex<float> temp[12][12]; 
  
    for (int i=0; i<12; i++) 
    { 
        for (int j=0; j<12; j++) 
        { 
            // Get cofactor of A[i][j] 
            getCofactor(A, temp, i, j, 12); 
  
            // sign of adj[j][i] positive if sum of row 
            // and column indexes is even. 
            sign = ((i+j)%2==0)? 1: -1; 
  
            // Interchanging rows and columns to get the 
            // transpose of the cofactor matrix 
            adj[j][i] = (sign)*(determinant(temp, 11)); 
        } 
    } 
} 

// Function to calculate and store inverse, returns false if 
// matrix is singular 
void inverse(complex<float> A[][12], complex<float> inverse[][12]) 
{ cout<<'a';
    // Find determinant of A[][] 
    complex<float> det = determinant(A, 12); 

    // Find adjoint 
    complex<float> adj[12][12]; 
    adjoint(A, adj); 
  
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
    for (int i=0; i<12; i++) 
        {for (int j=0; j<12; j++) 
            {inverse[i][j] = adj[i][j]/det;}} 
  
    return ;
} 





void matmult(vector<vector<complex<float> > >   &A , vector<vector<complex<float> > >   &B, vector<vector<complex<float> > >   &C, int Am, int An, int Bm, int Bn )
{
  if(An != Bm){
    cout<<"error";
        return;}
else
{
for(int m=0;m < Am; m++)
    {for(int n=0; n <Bn ; n++)
        {complex<float> count(0.0,0.0) ;
            for (int k = 0; k<An; k++)
                {count = count + A[m][k]*B[k][n];}
            C[m][n] = count;}}}
            return;
}

void transpose(vector<vector<complex<float> > >  &A, int m, int n, vector<vector<complex<float> > >  &B)
{for(int i =0 ; i < m ; i++)
for(int j = 0; j < n ;j++)
{B[j][i] = conj(A[i][j]);}
}

int main(){
    cout<<"start";
  vector<vector<vector<complex<float> > > > ksp(256,vector<vector<complex<float> > > (256,vector<complex<float> >(2)));

for(int i = 0; i< 256; i++)
{for(int j= 0; j< 256; j++)
{complex <float> a(real1[i][j],imag1[i][j]);
  ksp[i][j][0] = a;
complex <float> b(real2[i][j],imag2[i][j]);
  ksp[i][j][1] = b;
}
}

complex<float> P(0.0,0.0);
vector<vector<vector<complex<float> > > > weights(12,vector<vector<complex<float> > > (2,vector<complex<float> >(256)));
vector<vector<vector<complex<float> > > > kspace(256,vector<vector<complex<float> > > (256,vector<complex<float> >(2)));
int i = 0;
while (i < 100)
    {if (i % 2 ==0)
        {for(int j = 0; j < 256 ; j++)
            {//kspace[i][j][0] = ksp[i][j][0];
               // kspace[i][j][1] = ksp[i][j][1];
          kspace[i][j][0] = 0.0;
                          kspace[i][j][1] = 0.0;
            }
        }
    else
        {for(int j = 0; j < 256 ; j++)
            {kspace[i][j][0] = ksp[i][j][0];
                kspace[i][j][1] = ksp[i][j][1];}
            }
    i = i + 1;}
for (int i = 100; i< 156;i++)
    {for(int j = 0; j<256;j++)
        {for(int l = 0; l<2;l++)
            {kspace[i][j][l] = ksp[i][j][l];}}}
/*for j=101:156
    kspace(j,:,:) = ksp(j,:,:);*/


i = 156;
while (i < 256)
    {if (i % 2 ==1)
        {for(int j = 0; j < 256 ; j++)
            {
          //kspace[i][j][0] = ksp[i][j][0];
               // kspace[i][j][1] = ksp[i][j][1];
                kspace[i][j][0] = 0.0;
                                          kspace[i][j][1] = 0.0;}
        }
    else
        {for(int j = 0; j < 256 ; j++)
            {kspace[i][j][0] = ksp[i][j][0];
                kspace[i][j][1] = ksp[i][j][1];}
            }
    i = i + 1;}

vector<vector<complex<float> > >  input (56,vector<complex<float> >(12));
vector<vector<complex<float> > >  output(56,vector<complex<float> >(2));
for(int m=2;m<3;m++)
    { for(int j =100; j<156;j++){
             if ((m >0) && (m < 255))
             {output[j-100][0] = kspace[j][m][0];
                output[j-100][1] = kspace[j][m][1];
             input[j-100][0] = kspace[j-1][m-1][0];
             input[j-100][1] = kspace[j-1][m][0];
             input[j-100][2] = kspace[j-1][m+1][0];
             input[j-100][3] = kspace[j+1][m-1][0];
             input[j-100][4] = kspace[j+1][m][0];
             input[j-100][5] = kspace[j+1][m+1][0];
             input[j-100][6] = kspace[j-1][m-1][1];
             input[j-100][7]= kspace[j-1][m][1];
             input[j-100][8]= kspace[j-1][m+1][1];
             input[j-100][9] = kspace[j+1][m-1][1];
             input[j-100][10] = kspace[j+1][m][1];
             input[j-100][11] = kspace[j+1][m+1][1];}
             
        else if (m ==0)
                   {output[j-100][0] = kspace[j][m][0];
                output[j-100][1] = kspace[j][m][1];
             input[j-100][0] = P;
             input[j-100][1] = kspace[j-1][m][0];
             input[j-100][2] = kspace[j-1][m+1][0];
             input[j-100][3] = P;
             input[j-100][4] = kspace[j+1][m][0];
             input[j-100][5] = kspace[j+1][m+1][0];
             input[j-100][6] = P;
             input[j-100][7]= kspace[j-1][m][1];
             input[j-100][8]= kspace[j-1][m+1][1];
             input[j-100][9] = P;
             input[j-100][10] = kspace[j+1][m][1];
             input[j-100][11] = kspace[j+1][m+1][1];}
             
        else if (m==255)
                {output[j-100][0] = kspace[j][m][0];
                output[j-100][1] = kspace[j][m][1];
             input[j-100][0] = kspace[j-1][m-1][0];
             input[j-100][1] = kspace[j-1][m][0];
             input[j-100][2] = P;
             input[j-100][3] = kspace[j+1][m-1][0];
             input[j-100][4] = kspace[j+1][m][0];
             input[j-100][5] = P;
             input[j-100][6] = kspace[j-1][m-1][1];
             input[j-100][7]= kspace[j-1][m][1];
             input[j-100][8]= P;
             input[j-100][9] = kspace[j+1][m-1][1];
             input[j-100][10] = kspace[j+1][m][1];
             input[j-100][11] = P;}}
cout<<"c";
if (m  ==0 || m == 255)
{for(int a=0;a<12;a++)
{
    weights[a][0][m] = 0.0;
       weights[a][1][m] = 0.0;
}}
else
{vector<vector<complex<float> > >  conju(12,vector<complex<float> >(56));
transpose(input,56,12,conju);
vector<vector<complex<float> > >  mult(12,vector<complex<float> >(12));
matmult(conju,input,mult,12,56,56,12);
cout<<"b";
vector<vector<complex<float> > >  inv(12,vector<complex<float> >(12));
complex<float> inver[12][12];
complex<float> multi[12][12];
for(int i=0; i < 12; i++)
{for(int j=0; j< 12; j++)
{multi[i][j] = mult[i][j];
cout<<multi[i][j]<<" ";}
cout<<endl;}
/*cout<<m;
cout<<multi[8][5]<<endl;*/
/*cout<<multi[4][1]<<" ";
cout<<multi[8][8]<<" "<<m<<endl;*/
GJinverse(multi);
cout<<"a";
for(int i=0; i < 12; i++)
{for(int j=0; j< 12; j++)
{inv[i][j] = multi[i][j];
cout<<inv[i][j]<<" ";}
cout<<endl;}
/*cout<<multi[4][1]<<" ";
cout<<multi[8][8]<<" "<<m<<endl;*/
vector<vector<complex<float> > >  mult2(12,vector<complex<float> >(2));
matmult(conju,output,mult2,12,56,56,2);
vector<vector<complex<float> > >  fin(12,vector<complex<float> >(2));
matmult(inv,mult2,fin,12,12,12,2);
    //ans = (inv((input')*input))*((input')*(output)));
for(int a=0;a<12;a++)
{
    weights[a][0][m] = fin[a][0];
       weights[a][1][m] = fin[a][1];
}
cout<<"yes"<<endl;
/*cout<<weights[3][0][m]<<" ";
cout<<weights[8][1][m]<<" "<<m<<endl;*/

}

}
                                                      


vector<vector<complex<float> > >  in(1,vector<complex<float> >(12));
vector<vector<vector<complex<float> > > > final(256,vector<vector<complex<float> > > (256,vector<complex<float> >(2)));
for(int m=0; m <256;m++){
for ( int i = 0; i < 100 ; i++){
        if (i % 2 ==1)
            {for(int j=0; j< 256; j++)
                {for (int m = 0; m< 2; m++)
                    {final[i][j][m] = kspace [i][j][m];}}
            }
        else
              { if (i == 0)
              {in[0][1] = 0;
                in[0][7] = 0;
                in[0][4] = kspace[i+1][m][0];
                in[0][10] = kspace[i+1][m][1];
                in[0][0] = 0;
                in[0][9] = 0;
                in[0][6] = 0;
                in[0][3] = 0;
              in[0][2] = 0;
               in[0][5] = 0;
                in[0][8] = 0;
                in[0][11] = 0;
            } 
            else {
                in[0][1] = kspace[i-1][m][0];
                in[0][7]= kspace[i-1][m][1];
                in[0][4] = kspace[i+1][m][0];
                in[0][10] = kspace[i+1][m][1];

            if (m >0 && m < 255)
                {in[0][0] = kspace[i-1][m-1][0];
                in[0][9] = kspace[i+1][m-1][1];
                in[0][6] = kspace[i-1][m-1][1];
                in[0][3] = kspace[i+1][m-1][0];
               in[0][2] = kspace[i-1][m+1][0];
               in[0][5] = kspace[i+1][m+1][0];
                in[0][8] = kspace[i-1][m+1][1];
                in[0][11] = kspace[i+1][m+1][1];}
           else  if (m==0)
                {in[0][2] = kspace[i-1][m+1][0];
               in[0][5]= kspace[i+1][m+1][0];
                in[0][8] = kspace[i-1][m+1][1];
                in[0][11] = kspace[i+1][m+1][1];
                in[0][0] = 0;
                in[0][9] = 0;
                in[0][6]= 0;
                in[0][3] = 0;}
            else if (m == 255)
                {in[0][0] = kspace[i-1][m-1][0];
                in[0][9] = kspace[i+1][m-1][1];
                in[0][6] = kspace[i-1][m-1][1];
                in[0][3] = kspace[i+1][m-1][0];
                   in[0][2] = 0;
               in[0][5] = 0;
                in[0][8] = 0;
                in[0][11] = 0;}
            }
        }
    

for(int k=0;k<2;k++)
{complex<float> count(0.0,0.0);
    for(int p =0 ;p <12 ; p++)
{ count = count + (in[0][p]) * (weights[p][k][m]); }
      final[i][m][k] = count;
}
}

    for ( int i = 100; i<156;i++)
        {final[i][m][0] = kspace[i][m][0];
            final[i][m][1] = kspace[i][m][1];}


    for ( int i = 156; i<256;i++)
        {if (i % 2 ==0)
            {for(int j=0; j< 256; j++)
                {for (int z = 0; z< 2; z++)
                    {final[i][j][z] = kspace [i][j][z];}}
                }
        else
               {if (i == 255)
                {in[0][1] = kspace[i-1][m][0];
                in[0][7]= kspace[i-1][m][1];
                in[0][4] = 0;
                in[0][10] = 0;
                in[0][0] = 0;
                in[0][9] = 0;
                in[0][6] = 0;
                in[0][3] = 0;
               in[0][2] = 0;
               in[0][5] = 0;
                in[0][8] = 0;
                in[0][11] = 0;

                }
            else{
                in[0][1] = kspace[i-1][m][0];
                in[0][7]= kspace[i-1][m][1];
                in[0][4] = kspace[i+1][m][0];
                in[0][10] = kspace[i+1][m][1];
            if (m >0 && m < 255)
               {in[0][0] = kspace[i-1][m-1][0];
                in[0][9] = kspace[i+1][m-1][1];
                in[0][6] = kspace[i-1][m-1][1];
                in[0][3] = kspace[i+1][m-1][0];
               in[0][2] = kspace[i-1][m+1][0];
               in[0][5] = kspace[i+1][m+1][0];
                in[0][8] = kspace[i-1][m+1][1];
                in[0][11] = kspace[i+1][m+1][1];}
            else if (m==0)
                {in[0][2] = kspace[i-1][m+1][0];
               in[0][5]= kspace[i+1][m+1][0];
                in[0][8] = kspace[i-1][m+1][1];
                in[0][11] = kspace[i+1][m+1][1];
                in[0][0] = 0;
                in[0][9] = 0;
                in[0][6]= 0;
                in[0][3] = 0;}
            else if (m == 255)
                {in[0][0] = kspace[i-1][m-1][0];
                in[0][9] = kspace[i+1][m-1][1];
                in[0][6] = kspace[i-1][m-1][1];
                in[0][3] = kspace[i+1][m-1][0];
                   in[0][2] = 0;
               in[0][5] = 0;
                in[0][8] = 0;
                in[0][11] = 0;}}}



for(int k=0;k<2;k++)
{complex<float> count(0.0,0.0);
    for(int p =0 ;p <12 ; p++)
{ count = count + in[0][p]*weights[p][k][m]; }
      final[i][m][k] = count;
}
}
}

  //final(i,m,:) = (in)*(weights(:,:,m));
//combin = sqrt((pehla.*pehla + dusra.*dusra));
//final1 = (ifft((pehla)));
//final2 = (ifft((dusra)));
vector<vector<complex<float> > >  combine(256,vector<complex<float> >(256));
for(int i=0;i<256;i++)
{for(int j=0;j<256;j++)
{combine[i][j] = sqrt((final[i][j][0])*(final[i][j][0]) + (final[i][j][1])*(final[i][j][1])) ;
}
}

for(int j=0;j<256;j++)
{cout<<combine[5][j]<<endl;
}



cout<<"finish"<<endl;


}

/*#include<iostream>
using namespace std;
int main()
{
   int i,j,k,n;
   float a[100][200],t;
   n =12;
   for(i=0;i<n;i++)
   {
      for(j=n;j<2*n;j++)
      {
          if(i==j-n)
             a[i][j]=1;
         else
             a[i][j]=0;
       }
   }
   for(i=0;i<n;i++)
   {
      t=a[i][i];
      for(j=i;j<2*n;j++)
          a[i][j]=a[i][j]/t;
      for(j=0;j<n;j++)
      {
         if(i!=j)
         {
            t=a[j][i];
            for(k=0;k<2*n;k++)
                a[j][k]=a[j][k]-t*a[i][k];
          }
      }
   }
   cout<<"\n\nInverse matrix\n\n";
   for(i=0;i<n;i++)
   {
      for(j=n;j<2*n;j++)
         cout<<"\t"<<a[i][j];
      cout<<"\n";
    }
return 0;
}*/

