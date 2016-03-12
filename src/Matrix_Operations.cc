#include "Matrix_Operations.hh"

#include <iostream>
#include <vector>

#include "Check.hh"

using namespace std;

Matrix_Operations::
Matrix_Operations()
{
}

// dot product a . b
void Matrix_Operations::
dot(vector<double> const &a_data,
    vector<double> const &b_data,
    double &x_data,
    unsigned number_of_elements)
{
    Check(a_data.size() == number_of_elements, "a size");
    Check(b_data.size() == number_of_elements, "b size");
    
    double sum = 0;

    for (unsigned i=0; i<number_of_elements; ++i)
    {
        sum += a_data[i]*b_data[i];
    }

    x_data = sum;
}


// cross product A x B
void Matrix_Operations::
cross(vector<double> const &a_data,
      vector<double> const &b_data,
      vector<double> &x_data,
      unsigned number_of_elements)
{
    Check(number_of_elements == 3, "size must be 3 for cross product");
    Check(a_data.size() == number_of_elements, "a size");
    Check(b_data.size() == number_of_elements, "b size");
    Check(x_data.size() == number_of_elements, "x size");
    
    x_data[0] = a_data[1]*b_data[2] - a_data[2]*b_data[1];
    x_data[1] = a_data[2]*b_data[0] - a_data[0]*b_data[2];
    x_data[2] = a_data[0]*b_data[1] - a_data[1]*b_data[0];
}

// multiplies matrices A (n x m) and B (m x p)
void Matrix_Operations::
multiply(vector<double> const &a_data,
         vector<double> const &b_data,
         vector<double> &x_data,
         unsigned n,
         unsigned m,
         unsigned p)
{
    Check(a_data.size() == n * m, "a size");
    Check(b_data.size() == m * p, "b size");
    Check(x_data.size() == n * p, "x size");
    
    for (unsigned i=0; i<n; ++i)
    {
        for (unsigned j=0; j<p; ++j)
        {
            double sum = 0;
            
            for (unsigned k=0; k<m; ++k)
            {
                unsigned index1 = k + m*i;
                unsigned index2 = j + p*k;
                
                sum += a_data[index1]*b_data[index2];
            }
            
            unsigned index = j + p*i;

            x_data[index] = sum;
        }
    }
}

// solves Ax=B analytically for a 2x2 matrix
void Matrix_Operations::
solve(vector<double> const &a,
      vector<double> const &b,
      vector<double> &x,
      unsigned number_of_elements)
{
    Check(a.size() == number_of_elements*number_of_elements, "a size");
    Check(b.size() == number_of_elements, "b size");
    Check(x.size() == number_of_elements, "x size");
    Check(number_of_elements < 7, "max number_of_elements = 6");
    
    switch(number_of_elements)
    {
    case 1:
        x[0] = b[0] / a[0];
        return;
    case 2:
    {
        double o1 =  a[2] * a[1] - a[0] * a[3];
        x[0] = (a[1] * b[1] - a[3] * b[0]) / o1;
        x[1] = (a[2] * b[0] - a[0] * b[1]) / o1;
    }
    return;
    case 3:
    {
        double o5=a[6];
        double o3=a[2];
        double o8=a[5];
        double o11=a[7];
        double o7=a[1];
        double o10=a[3];
        double o13=a[0];
        double o4=a[4];
        double o15=a[8];
        double o20=b[0];
        double o23=b[1];
        double o26=b[2];
        double o6=o3*o4*o5;
        double o9=-(o5*o7*o8);
        double o12=-(o10*o11*o3);
        double o14=o11*o13*o8;
        double o16=o10*o15*o7;
        double o17=-(o13*o15*o4);
        double o18=o12+o14+o16+o17+o6+o9;
        double o19=1/o18;
        x[0]=o19*(o10*o15*o23-o10*o11*o26-o15*o20*o4+o26*o4*o5+o11*o20*o8-o23*o5*o8);
        x[1]=o19*(-(o13*o15*o23)+o11*o13*o26-o11*o20*o3+o23*o3*o5+o15*o20*o7-o26*o5*o7);
        x[2]=o19*(-(o10*o23*o3)-o13*o26*o4+o20*o3*o4+o10*o26*o7+o13*o23*o8-o20*o7*o8);
    }
    return;
    case 4:
    {
        double o5=a[9];
        double o6=a[12];
        double o3=a[3];
        double o9=a[7];
        double o12=a[10];
        double o8=a[2];
        double o11=a[5];
        double o14=a[1];
        double o4=a[6];
        double o16=a[11];
        double o19=a[8];
        double o20=a[13];
        double o23=a[4];
        double o25=a[0];
        double o29=a[14];
        double o36=a[15];
        double o45=b[0];
        double o52=b[1];
        double o59=b[2];
        double o66=b[3];
        double o7=-(o3*o4*o5*o6);
        double o10=o5*o6*o8*o9;
        double o13=o11*o12*o3*o6;
        double o15=-(o12*o14*o6*o9);
        double o17=-(o11*o16*o6*o8);
        double o18=o14*o16*o4*o6;
        double o21=o19*o20*o3*o4;
        double o22=-(o19*o20*o8*o9);
        double o24=-(o12*o20*o23*o3);
        double o26=o12*o20*o25*o9;
        double o27=o16*o20*o23*o8;
        double o28=-(o16*o20*o25*o4);
        double o30=-(o11*o19*o29*o3);
        double o31=o14*o19*o29*o9;
        double o32=o23*o29*o3*o5;
        double o33=-(o25*o29*o5*o9);
        double o34=-(o14*o16*o23*o29);
        double o35=o11*o16*o25*o29;
        double o37=o11*o19*o36*o8;
        double o38=-(o14*o19*o36*o4);
        double o39=-(o23*o36*o5*o8);
        double o40=o25*o36*o4*o5;
        double o41=o12*o14*o23*o36;
        double o42=-(o11*o12*o25*o36);
        double o43=o10+o13+o15+o17+o18+o21+o22+o24+o26+o27+o28+o30+o31+o32+o33+o34+o35+o37+o38+o39+o40+o41+o42+o7;
        double o44=1/o43;
        double o134=-(o14*o23);
        double o135=o11*o25;
        double o136=o134+o135;
        double o142=-(o23*o8);
        double o143=o25*o4;
        double o144=o142+o143;
        double o130=-(o14*o19);
        double o131=o25*o5;
        double o132=o130+o131;
        double o127=-(o23*o3);
        double o128=o25*o9;
        double o129=o127+o128;
        double o145=-(o14*o6);
        double o146=o20*o25;
        double o147=o145+o146;
        double o148=-(o144*o147);
        double o149=-(o6*o8);
        double o150=o25*o29;
        double o151=o149+o150;
        double o152=o136*o151;
        double o153=o148+o152;
        double o155=-(o132*o144);
        double o156=-(o19*o8);
        double o157=o12*o25;
        double o158=o156+o157;
        double o159=o136*o158;
        double o160=o155+o159;
        double o170=-(o23*o45);
        double o171=o25*o52;
        double o172=o170+o171;
        x[0]=o44*(o11*o16*o29*o45-o11*o12*o36*o45-o16*o20*o4*o45+o36*o4*o45*o5-o14*o16*o29*o52-o12*o20*o3*o52+o12*o14*o36*o52+o29*o3*o5*o52-o11*o29*o3*o59+o20*o3*o4*o59-o14*o36*o4*o59+o11*o12*o3*o66+o14*o16*o4*o66-o3*o4*o5*o66+o16*o20*o52*o8-o36*o5*o52*o8+o11*o36*o59*o8-o11*o16*o66*o8+o12*o20*o45*o9-o29*o45*o5*o9+o14*o29*o59*o9-o12*o14*o66*o9-o20*o59*o8*o9+o5*o66*o8*o9);
        x[1]=o44*(-(o16*o23*o29*o45)+o12*o23*o36*o45-o19*o36*o4*o45+o16*o25*o29*o52-o19*o29*o3*o52-o12*o25*o36*o52+o23*o29*o3*o59+o25*o36*o4*o59+o16*o4*o45*o6+o12*o3*o52*o6-o3*o4*o59*o6-o12*o23*o3*o66-o16*o25*o4*o66+o19*o3*o4*o66+o19*o36*o52*o8-o23*o36*o59*o8-o16*o52*o6*o8+o16*o23*o66*o8+o19*o29*o45*o9-o25*o29*o59*o9-o12*o45*o6*o9+o12*o25*o66*o9+o59*o6*o8*o9-o19*o66*o8*o9);
        x[2]=o44*(o16*o20*o23*o45+o11*o19*o36*o45-o23*o36*o45*o5-o16*o20*o25*o52+o19*o20*o3*o52-o14*o19*o36*o52+o25*o36*o5*o52-o20*o23*o3*o59+o14*o23*o36*o59-o11*o25*o36*o59-o11*o16*o45*o6+o14*o16*o52*o6-o3*o5*o52*o6+o11*o3*o59*o6-o14*o16*o23*o66+o11*o16*o25*o66-o11*o19*o3*o66+o23*o3*o5*o66-o19*o20*o45*o9+o20*o25*o59*o9+o45*o5*o6*o9-o14*o59*o6*o9+o14*o19*o66*o9-o25*o5*o66*o9);
        x[3]=(-(o153*(-(o132*o172)+o136*(-(o19*o45)+o25*o59)))+o160*(-(o147*o172)+o136*(-(o45*o6)+o25*o66)))/(-(o153*(-(o129*o132)+o136*(o16*o25-o19*o3)))+o160*(-(o129*o147)+o136*(o25*o36-o3*o6)));
    }
    return;
    case 5:
    {
        double o5=a[11];
        double o6=a[15];
        double o3=a[3];
        double o9=a[8];
        double o12=a[12];
        double o8=a[2];
        double o11=a[6];
        double o14=a[1];
        double o4=a[7];
        double o16=a[13];
        double o19=a[10];
        double o20=a[16];
        double o23=a[5];
        double o25=a[0];
        double o29=a[17];
        double o36=a[18];
        double o45=b[0];
        double o52=b[1];
        double o59=b[2];
        double o66=b[3];
        double o7=-(o3*o4*o5*o6);
        double o10=o5*o6*o8*o9;
        double o13=o11*o12*o3*o6;
        double o15=-(o12*o14*o6*o9);
        double o17=-(o11*o16*o6*o8);
        double o18=o14*o16*o4*o6;
        double o21=o19*o20*o3*o4;
        double o22=-(o19*o20*o8*o9);
        double o24=-(o12*o20*o23*o3);
        double o26=o12*o20*o25*o9;
        double o27=o16*o20*o23*o8;
        double o28=-(o16*o20*o25*o4);
        double o30=-(o11*o19*o29*o3);
        double o31=o14*o19*o29*o9;
        double o32=o23*o29*o3*o5;
        double o33=-(o25*o29*o5*o9);
        double o34=-(o14*o16*o23*o29);
        double o35=o11*o16*o25*o29;
        double o37=o11*o19*o36*o8;
        double o38=-(o14*o19*o36*o4);
        double o39=-(o23*o36*o5*o8);
        double o40=o25*o36*o4*o5;
        double o41=o12*o14*o23*o36;
        double o42=-(o11*o12*o25*o36);
        double o43=o10+o13+o15+o17+o18+o21+o22+o24+o26+o27+o28+o30+o31+o32+o33+o34+o35+o37+o38+o39+o40+o41+o42+o7;
        double o44=1/o43;
        double o75=a[4];
        double o77=a[9];
        double o81=a[14];
        double o96=a[19];
        double o111=-(o14*o23);
        double o112=o11*o25;
        double o113=o111+o112;
        double o119=-(o23*o8);
        double o120=o25*o4;
        double o121=o119+o120;
        double o107=-(o14*o19);
        double o108=o25*o5;
        double o109=o107+o108;
        double o104=-(o23*o75);
        double o105=o25*o77;
        double o106=o104+o105;
        double o122=-(o14*o6);
        double o123=o20*o25;
        double o124=o122+o123;
        double o155=a[20];
        double o132=-(o109*o121);
        double o133=-(o19*o8);
        double o134=o12*o25;
        double o135=o133+o134;
        double o136=o113*o135;
        double o137=o132+o136;
        double o146=-(o23*o3);
        double o147=o25*o9;
        double o148=o146+o147;
        double o156=-(o14*o155);
        double o157=a[21];
        double o158=o157*o25;
        double o159=o156+o158;
        double o149=-(o109*o148);
        double o150=-(o19*o3);
        double o151=o16*o25;
        double o152=o150+o151;
        double o153=o113*o152;
        double o154=o149+o153;
        double o125=-(o121*o124);
        double o126=-(o6*o8);
        double o127=o25*o29;
        double o128=o126+o127;
        double o129=o113*o128;
        double o130=o125+o129;
        double o110=-(o106*o109);
        double o114=-(o19*o75);
        double o115=o25*o81;
        double o116=o114+o115;
        double o117=o113*o116;
        double o118=o110+o117;
        double o160=-(o121*o159);
        double o161=-(o155*o8);
        double o162=a[22];
        double o163=o162*o25;
        double o164=o161+o163;
        double o165=o113*o164;
        double o166=o160+o165;
        double o167=-(o154*o166);
        double o168=-(o148*o159);
        double o169=-(o155*o3);
        double o170=a[23];
        double o171=o170*o25;
        double o172=o169+o171;
        double o173=o113*o172;
        double o174=o168+o173;
        double o175=o137*o174;
        double o176=o167+o175;
        double o200=-(o23*o45);
        double o201=o25*o52;
        double o202=o200+o201;
        double o178=-(o130*o154);
        double o179=-(o124*o148);
        double o180=-(o3*o6);
        double o181=o25*o36;
        double o182=o180+o181;
        double o183=o113*o182;
        double o184=o179+o183;
        double o185=o137*o184;
        double o186=o178+o185;
        double o203=-(o109*o202);
        double o204=-(o19*o45);
        double o205=o25*o59;
        double o206=o204+o205;
        double o207=o113*o206;
        double o208=o203+o207;
        double o131=-(o118*o130);
        double o138=-(o106*o124);
        double o139=-(o6*o75);
        double o140=o25*o96;
        double o141=o139+o140;
        double o142=o113*o141;
        double o143=o138+o142;
        double o144=o137*o143;
        double o145=o131+o144;
        double o177=-(o145*o176);
        double o187=-(o118*o166);
        double o188=-(o106*o159);
        double o189=-(o155*o75);
        double o190=a[24];
        double o191=o190*o25;
        double o192=o189+o191;
        double o193=o113*o192;
        double o194=o188+o193;
        double o195=o137*o194;
        double o196=o187+o195;
        double o197=o186*o196;
        double o198=o177+o197;
        double o199=1/o198;
        double o209=-(o130*o208);
        double o210=-(o124*o202);
        double o211=-(o45*o6);
        double o212=o25*o66;
        double o213=o211+o212;
        double o214=o113*o213;
        double o215=o210+o214;
        double o216=o137*o215;
        double o217=o209+o216;
        double o218=-(o176*o217);
        double o219=-(o166*o208);
        double o220=-(o159*o202);
        double o221=-(o155*o45);
        double o222=b[4];
        double o223=o222*o25;
        double o224=o221+o223;
        double o225=o113*o224;
        double o226=o220+o225;
        double o227=o137*o226;
        double o228=o219+o227;
        double o229=o186*o228;
        double o230=o218+o229;
        double o339=1/o186;
        x[0]=o44*(o11*o16*o29*o45-o11*o12*o36*o45-o16*o20*o4*o45+o36*o4*o45*o5-o14*o16*o29*o52-o12*o20*o3*o52+o12*o14*o36*o52+o29*o3*o5*o52-o11*o29*o3*o59+o20*o3*o4*o59-o14*o36*o4*o59+o11*o12*o3*o66+o14*o16*o4*o66-o3*o4*o5*o66+o16*o20*o52*o8-o36*o5*o52*o8+o11*o36*o59*o8-o11*o16*o66*o8+o12*o20*o45*o9-o29*o45*o5*o9+o14*o29*o59*o9-o12*o14*o66*o9-o20*o59*o8*o9+o5*o66*o8*o9)-o199*o230*o44*(o11*o16*o29*o75-o11*o12*o36*o75-o16*o20*o4*o75+o36*o4*o5*o75-o14*o16*o29*o77-o12*o20*o3*o77+o12*o14*o36*o77+o29*o3*o5*o77+o16*o20*o77*o8-o36*o5*o77*o8-o11*o29*o3*o81+o20*o3*o4*o81-o14*o36*o4*o81+o11*o36*o8*o81+o12*o20*o75*o9-o29*o5*o75*o9+o14*o29*o81*o9-o20*o8*o81*o9+o11*o12*o3*o96+o14*o16*o4*o96-o3*o4*o5*o96-o11*o16*o8*o96-o12*o14*o9*o96+o5*o8*o9*o96);
        x[1]=o44*(-(o16*o23*o29*o45)+o12*o23*o36*o45-o19*o36*o4*o45+o16*o25*o29*o52-o19*o29*o3*o52-o12*o25*o36*o52+o23*o29*o3*o59+o25*o36*o4*o59+o16*o4*o45*o6+o12*o3*o52*o6-o3*o4*o59*o6-o12*o23*o3*o66-o16*o25*o4*o66+o19*o3*o4*o66+o19*o36*o52*o8-o23*o36*o59*o8-o16*o52*o6*o8+o16*o23*o66*o8+o19*o29*o45*o9-o25*o29*o59*o9-o12*o45*o6*o9+o12*o25*o66*o9+o59*o6*o8*o9-o19*o66*o8*o9)-o199*o230*o44*(-(o16*o23*o29*o75)+o12*o23*o36*o75-o19*o36*o4*o75+o16*o4*o6*o75+o16*o25*o29*o77-o19*o29*o3*o77-o12*o25*o36*o77+o12*o3*o6*o77+o19*o36*o77*o8-o16*o6*o77*o8+o23*o29*o3*o81+o25*o36*o4*o81-o3*o4*o6*o81-o23*o36*o8*o81+o19*o29*o75*o9-o12*o6*o75*o9-o25*o29*o81*o9+o6*o8*o81*o9-o12*o23*o3*o96-o16*o25*o4*o96+o19*o3*o4*o96+o16*o23*o8*o96+o12*o25*o9*o96-o19*o8*o9*o96);
        x[2]=o44*(o16*o20*o23*o45+o11*o19*o36*o45-o23*o36*o45*o5-o16*o20*o25*o52+o19*o20*o3*o52-o14*o19*o36*o52+o25*o36*o5*o52-o20*o23*o3*o59+o14*o23*o36*o59-o11*o25*o36*o59-o11*o16*o45*o6+o14*o16*o52*o6-o3*o5*o52*o6+o11*o3*o59*o6-o14*o16*o23*o66+o11*o16*o25*o66-o11*o19*o3*o66+o23*o3*o5*o66-o19*o20*o45*o9+o20*o25*o59*o9+o45*o5*o6*o9-o14*o59*o6*o9+o14*o19*o66*o9-o25*o5*o66*o9)-o199*o230*o44*(o16*o20*o23*o75+o11*o19*o36*o75-o23*o36*o5*o75-o11*o16*o6*o75-o16*o20*o25*o77+o19*o20*o3*o77-o14*o19*o36*o77+o25*o36*o5*o77+o14*o16*o6*o77-o3*o5*o6*o77-o20*o23*o3*o81+o14*o23*o36*o81-o11*o25*o36*o81+o11*o3*o6*o81-o19*o20*o75*o9+o5*o6*o75*o9+o20*o25*o81*o9-o14*o6*o81*o9-o14*o16*o23*o96+o11*o16*o25*o96-o11*o19*o3*o96+o23*o3*o5*o96+o14*o19*o9*o96-o25*o5*o9*o96);
        x[3]=o217*o339-o145*o199*o230*o339;
        x[4]=o199*o230;
    }
    break;
    case 6:
    {
        double o5=a[13];
        double o6=a[18];
        double o3=a[3];
        double o9=a[9];
        double o12=a[14];
        double o8=a[2];
        double o11=a[7];
        double o14=a[1];
        double o4=a[8];
        double o16=a[15];
        double o19=a[12];
        double o20=a[19];
        double o23=a[6];
        double o25=a[0];
        double o29=a[20];
        double o36=a[21];
        double o45=b[0];
        double o52=b[1];
        double o59=b[2];
        double o66=b[3];
        double o7=-(o3*o4*o5*o6);
        double o10=o5*o6*o8*o9;
        double o13=o11*o12*o3*o6;
        double o15=-(o12*o14*o6*o9);
        double o17=-(o11*o16*o6*o8);
        double o18=o14*o16*o4*o6;
        double o21=o19*o20*o3*o4;
        double o22=-(o19*o20*o8*o9);
        double o24=-(o12*o20*o23*o3);
        double o26=o12*o20*o25*o9;
        double o27=o16*o20*o23*o8;
        double o28=-(o16*o20*o25*o4);
        double o30=-(o11*o19*o29*o3);
        double o31=o14*o19*o29*o9;
        double o32=o23*o29*o3*o5;
        double o33=-(o25*o29*o5*o9);
        double o34=-(o14*o16*o23*o29);
        double o35=o11*o16*o25*o29;
        double o37=o11*o19*o36*o8;
        double o38=-(o14*o19*o36*o4);
        double o39=-(o23*o36*o5*o8);
        double o40=o25*o36*o4*o5;
        double o41=o12*o14*o23*o36;
        double o42=-(o11*o12*o25*o36);
        double o43=o10+o13+o15+o17+o18+o21+o22+o24+o26+o27+o28+o30+o31+o32+o33+o34+o35+o37+o38+o39+o40+o41+o42+o7;
        double o44=1/o43;
        double o75=a[4];
        double o77=a[10];
        double o81=a[16];
        double o96=a[22];
        double o111=-(o14*o23);
        double o112=o11*o25;
        double o113=o111+o112;
        double o119=-(o23*o8);
        double o120=o25*o4;
        double o121=o119+o120;
        double o107=-(o14*o19);
        double o108=o25*o5;
        double o109=o107+o108;
        double o104=-(o23*o75);
        double o105=o25*o77;
        double o106=o104+o105;
        double o122=-(o14*o6);
        double o123=o20*o25;
        double o124=o122+o123;
        double o155=a[24];
        double o132=-(o109*o121);
        double o133=-(o19*o8);
        double o134=o12*o25;
        double o135=o133+o134;
        double o136=o113*o135;
        double o137=o132+o136;
        double o146=-(o23*o3);
        double o147=o25*o9;
        double o148=o146+o147;
        double o156=-(o14*o155);
        double o157=a[25];
        double o158=o157*o25;
        double o159=o156+o158;
        double o149=-(o109*o148);
        double o150=-(o19*o3);
        double o151=o16*o25;
        double o152=o150+o151;
        double o153=o113*o152;
        double o154=o149+o153;
        double o125=-(o121*o124);
        double o126=-(o6*o8);
        double o127=o25*o29;
        double o128=o126+o127;
        double o129=o113*o128;
        double o130=o125+o129;
        double o110=-(o106*o109);
        double o114=-(o19*o75);
        double o115=o25*o81;
        double o116=o114+o115;
        double o117=o113*o116;
        double o118=o110+o117;
        double o160=-(o121*o159);
        double o161=-(o155*o8);
        double o162=a[26];
        double o163=o162*o25;
        double o164=o161+o163;
        double o165=o113*o164;
        double o166=o160+o165;
        double o167=-(o154*o166);
        double o168=-(o148*o159);
        double o169=-(o155*o3);
        double o170=a[27];
        double o171=o170*o25;
        double o172=o169+o171;
        double o173=o113*o172;
        double o174=o168+o173;
        double o175=o137*o174;
        double o176=o167+o175;
        double o200=-(o23*o45);
        double o201=o25*o52;
        double o202=o200+o201;
        double o178=-(o130*o154);
        double o179=-(o124*o148);
        double o180=-(o3*o6);
        double o181=o25*o36;
        double o182=o180+o181;
        double o183=o113*o182;
        double o184=o179+o183;
        double o185=o137*o184;
        double o186=o178+o185;
        double o203=-(o109*o202);
        double o204=-(o19*o45);
        double o205=o25*o59;
        double o206=o204+o205;
        double o207=o113*o206;
        double o208=o203+o207;
        double o232=a[5];
        double o234=a[11];
        double o238=a[17];
        double o253=a[23];
        double o76=o12*o20*o75*o9;
        double o78=-(o12*o20*o3*o77);
        double o79=-(o16*o20*o4*o75);
        double o80=o16*o20*o77*o8;
        double o82=o20*o3*o4*o81;
        double o83=-(o20*o8*o81*o9);
        double o84=-(o29*o5*o75*o9);
        double o85=o29*o3*o5*o77;
        double o86=o11*o16*o29*o75;
        double o87=-(o14*o16*o29*o77);
        double o88=-(o11*o29*o3*o81);
        double o89=o14*o29*o81*o9;
        double o90=o36*o4*o5*o75;
        double o91=-(o36*o5*o77*o8);
        double o92=-(o11*o12*o36*o75);
        double o93=o12*o14*o36*o77;
        double o94=o11*o36*o8*o81;
        double o95=-(o14*o36*o4*o81);
        double o97=-(o3*o4*o5*o96);
        double o98=o5*o8*o9*o96;
        double o99=o11*o12*o3*o96;
        double o100=-(o12*o14*o9*o96);
        double o101=-(o11*o16*o8*o96);
        double o102=o14*o16*o4*o96;
        double o103=o100+o101+o102+o76+o78+o79+o80+o82+o83+o84+o85+o86+o87+o88+o89+o90+o91+o92+o93+o94+o95+o97+o98+o99;
        double o131=-(o118*o130);
        double o138=-(o106*o124);
        double o139=-(o6*o75);
        double o140=o25*o96;
        double o141=o139+o140;
        double o142=o113*o141;
        double o143=o138+o142;
        double o144=o137*o143;
        double o145=o131+o144;
        double o177=-(o145*o176);
        double o187=-(o118*o166);
        double o188=-(o106*o159);
        double o189=-(o155*o75);
        double o190=a[28];
        double o191=o190*o25;
        double o192=o189+o191;
        double o193=o113*o192;
        double o194=o188+o193;
        double o195=o137*o194;
        double o196=o187+o195;
        double o197=o186*o196;
        double o198=o177+o197;
        double o199=1/o198;
        double o262=-(o23*o232);
        double o263=o234*o25;
        double o264=o262+o263;
        double o265=-(o109*o264);
        double o266=-(o19*o232);
        double o267=o238*o25;
        double o268=o266+o267;
        double o269=o113*o268;
        double o270=o265+o269;
        double o271=-(o130*o270);
        double o272=-(o124*o264);
        double o273=-(o232*o6);
        double o274=o25*o253;
        double o275=o273+o274;
        double o276=o113*o275;
        double o277=o272+o276;
        double o278=o137*o277;
        double o279=o271+o278;
        double o280=-(o176*o279);
        double o281=-(o166*o270);
        double o282=-(o159*o264);
        double o283=-(o155*o232);
        double o284=a[29];
        double o285=o25*o284;
        double o286=o283+o285;
        double o287=o113*o286;
        double o288=o282+o287;
        double o289=o137*o288;
        double o290=o281+o289;
        double o291=o186*o290;
        double o292=o280+o291;
        double o295=a[30];
        double o296=-(o14*o295);
        double o297=a[31];
        double o298=o25*o297;
        double o299=o296+o298;
        double o300=-(o121*o299);
        double o301=-(o295*o8);
        double o302=a[32];
        double o303=o25*o302;
        double o304=o301+o303;
        double o305=o113*o304;
        double o306=o300+o305;
        double o307=-(o154*o306);
        double o308=-(o148*o299);
        double o309=-(o295*o3);
        double o310=a[33];
        double o311=o25*o310;
        double o312=o309+o311;
        double o313=o113*o312;
        double o314=o308+o313;
        double o315=o137*o314;
        double o316=o307+o315;
        double o317=-(o145*o316);
        double o318=-(o118*o306);
        double o319=-(o106*o299);
        double o320=-(o295*o75);
        double o321=a[34];
        double o322=o25*o321;
        double o323=o320+o322;
        double o324=o113*o323;
        double o325=o319+o324;
        double o326=o137*o325;
        double o327=o318+o326;
        double o328=o186*o327;
        double o329=o317+o328;
        double o209=-(o130*o208);
        double o210=-(o124*o202);
        double o211=-(o45*o6);
        double o212=o25*o66;
        double o213=o211+o212;
        double o214=o113*o213;
        double o215=o210+o214;
        double o216=o137*o215;
        double o217=o209+o216;
        double o218=-(o176*o217);
        double o219=-(o166*o208);
        double o220=-(o159*o202);
        double o221=-(o155*o45);
        double o222=b[4];
        double o223=o222*o25;
        double o224=o221+o223;
        double o225=o113*o224;
        double o226=o220+o225;
        double o227=o137*o226;
        double o228=o219+o227;
        double o229=o186*o228;
        double o230=o218+o229;
        double o391=-(o12*o6*o75*o9);
        double o392=o12*o3*o6*o77;
        double o393=o16*o4*o6*o75;
        double o394=-(o16*o6*o77*o8);
        double o395=-(o3*o4*o6*o81);
        double o396=o6*o8*o81*o9;
        double o397=o19*o29*o75*o9;
        double o398=-(o19*o29*o3*o77);
        double o399=-(o16*o23*o29*o75);
        double o400=o16*o25*o29*o77;
        double o401=o23*o29*o3*o81;
        double o402=-(o25*o29*o81*o9);
        double o403=-(o19*o36*o4*o75);
        double o404=o19*o36*o77*o8;
        double o405=o12*o23*o36*o75;
        double o406=-(o12*o25*o36*o77);
        double o407=-(o23*o36*o8*o81);
        double o408=o25*o36*o4*o81;
        double o409=o19*o3*o4*o96;
        double o410=-(o19*o8*o9*o96);
        double o411=-(o12*o23*o3*o96);
        double o412=o12*o25*o9*o96;
        double o413=o16*o23*o8*o96;
        double o414=-(o16*o25*o4*o96);
        double o415=o391+o392+o393+o394+o395+o396+o397+o398+o399+o400+o401+o402+o403+o404+o405+o406+o407+o408+o409+o410+o411+o412+o413+o414;
        double o330=-(o292*o329);
        double o331=-(o279*o316);
        double o332=-(o270*o306);
        double o333=-(o264*o299);
        double o334=-(o232*o295);
        double o335=a[35];
        double o336=o25*o335;
        double o337=o334+o336;
        double o338=o113*o337;
        double o339=o333+o338;
        double o340=o137*o339;
        double o341=o332+o340;
        double o342=o186*o341;
        double o343=o331+o342;
        double o344=o198*o343;
        double o345=o330+o344;
        double o346=1/o345;
        double o347=-(o230*o329);
        double o348=-(o217*o316);
        double o349=-(o208*o306);
        double o350=-(o202*o299);
        double o351=-(o295*o45);
        double o352=b[5];
        double o353=o25*o352;
        double o354=o351+o353;
        double o355=o113*o354;
        double o356=o350+o355;
        double o357=o137*o356;
        double o358=o349+o357;
        double o359=o186*o358;
        double o360=o348+o359;
        double o361=o198*o360;
        double o362=o347+o361;
        double o473=o5*o6*o75*o9;
        double o474=-(o3*o5*o6*o77);
        double o475=-(o11*o16*o6*o75);
        double o476=o14*o16*o6*o77;
        double o477=o11*o3*o6*o81;
        double o478=-(o14*o6*o81*o9);
        double o479=-(o19*o20*o75*o9);
        double o480=o19*o20*o3*o77;
        double o481=o16*o20*o23*o75;
        double o482=-(o16*o20*o25*o77);
        double o483=-(o20*o23*o3*o81);
        double o484=o20*o25*o81*o9;
        double o485=o11*o19*o36*o75;
        double o486=-(o14*o19*o36*o77);
        double o487=-(o23*o36*o5*o75);
        double o488=o25*o36*o5*o77;
        double o489=o14*o23*o36*o81;
        double o490=-(o11*o25*o36*o81);
        double o491=-(o11*o19*o3*o96);
        double o492=o14*o19*o9*o96;
        double o493=o23*o3*o5*o96;
        double o494=-(o25*o5*o9*o96);
        double o495=-(o14*o16*o23*o96);
        double o496=o11*o16*o25*o96;
        double o497=o473+o474+o475+o476+o477+o478+o479+o480+o481+o482+o483+o484+o485+o486+o487+o488+o489+o490+o491+o492+o493+o494+o495+o496;
        double o529=1/o186;
        x[0]=-(o103*o199*o230*o44)+o44*(o11*o16*o29*o45-o11*o12*o36*o45-o16*o20*o4*o45+o36*o4*o45*o5-o14*o16*o29*o52-o12*o20*o3*o52+o12*o14*o36*o52+o29*o3*o5*o52-o11*o29*o3*o59+o20*o3*o4*o59-o14*o36*o4*o59+o11*o12*o3*o66+o14*o16*o4*o66-o3*o4*o5*o66+o16*o20*o52*o8-o36*o5*o52*o8+o11*o36*o59*o8-o11*o16*o66*o8+o12*o20*o45*o9-o29*o45*o5*o9+o14*o29*o59*o9-o12*o14*o66*o9-o20*o59*o8*o9+o5*o66*o8*o9)-o346*o362*(-(o103*o199*o292*o44)+o44*(o11*o16*o232*o29-o14*o16*o234*o29-o12*o20*o234*o3+o11*o12*o253*o3-o11*o238*o29*o3-o11*o12*o232*o36+o12*o14*o234*o36-o16*o20*o232*o4+o14*o16*o253*o4+o20*o238*o3*o4-o14*o238*o36*o4+o234*o29*o3*o5-o253*o3*o4*o5+o232*o36*o4*o5+o16*o20*o234*o8-o11*o16*o253*o8+o11*o238*o36*o8-o234*o36*o5*o8+o12*o20*o232*o9-o12*o14*o253*o9+o14*o238*o29*o9-o232*o29*o5*o9-o20*o238*o8*o9+o253*o5*o8*o9));
        x[1]=-(o199*o230*o415*o44)+o44*(-(o16*o23*o29*o45)+o12*o23*o36*o45-o19*o36*o4*o45+o16*o25*o29*o52-o19*o29*o3*o52-o12*o25*o36*o52+o23*o29*o3*o59+o25*o36*o4*o59+o16*o4*o45*o6+o12*o3*o52*o6-o3*o4*o59*o6-o12*o23*o3*o66-o16*o25*o4*o66+o19*o3*o4*o66+o19*o36*o52*o8-o23*o36*o59*o8-o16*o52*o6*o8+o16*o23*o66*o8+o19*o29*o45*o9-o25*o29*o59*o9-o12*o45*o6*o9+o12*o25*o66*o9+o59*o6*o8*o9-o19*o66*o8*o9)-o346*o362*(-(o199*o292*o415*o44)+o44*(-(o16*o23*o232*o29)+o16*o234*o25*o29-o12*o23*o253*o3-o19*o234*o29*o3+o23*o238*o29*o3+o12*o23*o232*o36-o12*o234*o25*o36-o16*o25*o253*o4+o19*o253*o3*o4-o19*o232*o36*o4+o238*o25*o36*o4+o12*o234*o3*o6+o16*o232*o4*o6-o238*o3*o4*o6+o16*o23*o253*o8+o19*o234*o36*o8-o23*o238*o36*o8-o16*o234*o6*o8+o12*o25*o253*o9+o19*o232*o29*o9-o238*o25*o29*o9-o12*o232*o6*o9-o19*o253*o8*o9+o238*o6*o8*o9));
        x[2]=-(o199*o230*o44*o497)+o44*(o16*o20*o23*o45+o11*o19*o36*o45-o23*o36*o45*o5-o16*o20*o25*o52+o19*o20*o3*o52-o14*o19*o36*o52+o25*o36*o5*o52-o20*o23*o3*o59+o14*o23*o36*o59-o11*o25*o36*o59-o11*o16*o45*o6+o14*o16*o52*o6-o3*o5*o52*o6+o11*o3*o59*o6-o14*o16*o23*o66+o11*o16*o25*o66-o11*o19*o3*o66+o23*o3*o5*o66-o19*o20*o45*o9+o20*o25*o59*o9+o45*o5*o6*o9-o14*o59*o6*o9+o14*o19*o66*o9-o25*o5*o66*o9)-o346*o362*(-(o199*o292*o44*o497)+o44*(o16*o20*o23*o232-o16*o20*o234*o25-o14*o16*o23*o253+o11*o16*o25*o253+o19*o20*o234*o3-o20*o23*o238*o3-o11*o19*o253*o3+o11*o19*o232*o36-o14*o19*o234*o36+o14*o23*o238*o36-o11*o238*o25*o36+o23*o253*o3*o5-o23*o232*o36*o5+o234*o25*o36*o5-o11*o16*o232*o6+o14*o16*o234*o6+o11*o238*o3*o6-o234*o3*o5*o6-o19*o20*o232*o9+o20*o238*o25*o9+o14*o19*o253*o9-o25*o253*o5*o9-o14*o238*o6*o9+o232*o5*o6*o9));
        x[3]=o217*o529-o145*o199*o230*o529-o346*o362*(o279*o529-o145*o199*o292*o529);
        x[4]=o199*o230-o199*o292*o346*o362;
        x[5]=o346*o362;
    }
    break;
    }
}
