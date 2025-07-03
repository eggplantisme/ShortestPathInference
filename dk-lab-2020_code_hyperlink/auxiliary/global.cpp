#define PI 3.14159265

#include <vector>
#include <cmath>
#include <assert.h>
#include <boost/math/special_functions/gamma.hpp>
#include "global.h"

int roundx(double x) // changed from round because causes conflict when compiling on unix

{
   if(x>0) {return (int)(2*x+1)/2;}
   else {return (int)(2*x-1)/2;}
}

double finite_alpha(double kp_max, double gamma)
{
    assert(gamma > 2);
    assert(kp_max > 0);
    return 1 -  std::exp( (2.-gamma) * std::log(kp_max) ) ;
}

void quickSort(std::vector<int>& numbers, std::vector<int>& addition)
{
//  cout << "int sort" << endl;
    int array_size = numbers.size();
    q_sort(numbers, addition, 0, array_size - 1);
    reverse(numbers);
    reverse(addition);
}

//**********************************************************************

#define IA 16807

#define IM 2147483647

#define AM (1.0/IM)

#define IQ 127773

#define IR 2836

#define NTAB 32

#define NDIV (1+(IM-1)/NTAB)

#define EPS 1.2e-7

#define RNMX (1.0-EPS)

float ran1(long *idum)


/*Minimal" random number generator of Park and Miller with Bays-Durham shuffle and added

safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint

values). Call with idum a negative integer to initialize; thereafter, do not alter idum between

successive deviates in a sequence. RNMX should approximate the largest floating value that is

less than 1.*/

{

    int j;

    long k;

    static long iy=0;

    static long iv[NTAB];

    float temp;

    if (*idum <= 0 || !iy) 

    {   /*Initialize.*/

        if (-(*idum) < 1) 

            *idum=1; /*Be sure to prevent idum = 0.*/

        else 

            *idum = -(*idum);

        for (j=NTAB+7;j>=0;j--) 

        {   /*Load the shuffle table (after 8 warm-ups).*/

            k=(*idum)/IQ;

            *idum=IA*(*idum-k*IQ)-IR*k;

            if (*idum < 0) *idum += IM;

            if (j < NTAB) iv[j] = *idum;

        }

        iy=iv[0];

    }

    k=(*idum)/IQ;   /*Start here when not initializing.*/

    *idum=IA*(*idum-k*IQ)-IR*k; /*Compute idum=(IA*idum) % IM without overlows by Schrage's method.*/

    if (*idum < 0)

        *idum += IM;

    j=iy/NDIV;  /*Will be in the range 0..NTAB-1.*/

    iy=iv[j];   /*Output previously stored value and refill the shuffle table.*/

    iv[j] = *idum;

    if ((temp=AM*iy) > RNMX) 

        return RNMX;    /* Because users don't expect endpoint values.*/

    else 

    return temp;

}

void shuffle_vector(std::vector<int>& v1, long *idum)
{
    int i, ctr;
    std::vector<int> v2;
    v2.resize(v1.size());
    int size = v1.size();
    
    for(i = 0; i < v2.size(); i++)
    {
        ctr = roundx(ran1(idum)*size);
        if (ctr == size)
        {
            ctr--;
        }
        v2[i] = v1[ctr];
        v1[ctr] = v1[size - 1];
        size--;
    }
    v1 = v2;
}

double distH2(double r1, double angle1, double r2, double angle2, double zeta)
{
    double angle = distS1(angle1, angle2);
    if (angle < 0.0000001)
    {
        //      cout << "1" << endl;
                return std::fabs(r1 - r2);
    }
    double dist;
    dist = (std::cosh(zeta*r1)*std::cosh(zeta*r2)) - (std::sinh(zeta*r1)*std::sinh(zeta*r2) * std::cos(angle));
    dist = std::acosh(dist);
    dist /= zeta;
    return dist;
    
}

double distH2approx(double r1, double angle1, double r2, double angle2, double zeta)
{
    double dist, angle = distS1(angle1, angle2);
    if (angle < 0.0000001)
    {
        //        cout << "1" << endl;
        return fabs(r1 - r2);
    }
    if (angle > 2. * sqrt(  exp(-zeta * r1) + exp(-zeta * r2) ))
    {
        return r1 + r2  + (2./zeta) * log (angle/2.);
    }
        dist = (cosh(zeta*r1)*cosh(zeta*r2)) - (sinh(zeta*r1)*sinh(zeta*r2) * cos(angle));
        return acosh(dist)/zeta;
}

int get_string_number(int num, char str[]) //gets the number which is the  numth in the str string
{
    int ctr1, ctr2;
    ctr1 = get_word_beginning(num, str);
    ctr2 = get_word_ending(ctr1, str);
    return get_number(ctr1, ctr2, str);
}

double distS1(double x1, double x2)
{
    return fabs(PI- fabs(PI - fabs(x1 - x2)));
}

int get_word_beginning(int num, char str[]) // finds the beginning of the num'th word
{
    int ctr1 = -1;
    int i;
    do
    {
        ctr1++;
    }
    while ( (str[ctr1] == ' ') || (str[ctr1] == '\t') || (str[ctr1] == '|') || (str[ctr1] == ',') || (str[ctr1] == '#') || (str[ctr1] == '<') || (str[ctr1] == '>') || (str[ctr1] == ')') || (str[ctr1] == '(') || (str[ctr1] == '"')  );
    
//  cout << "the beginning of the word " << 1 <<" : " << ctr1 << endl;
    for(i = 0; i < num - 1; i++)
    {
        //find the end of the word word
        do
        {
            ctr1++;
        }
        while( (str[ctr1] != ' ') && (str[ctr1] != '\t') && (str[ctr1] != '|') && (str[ctr1] != ',') && (str[ctr1] != '#') && (str[ctr1] != '>') && (str[ctr1] != '<') && (str[ctr1] != ')') && (str[ctr1] != '(') && (str[ctr1] != '"'));
//      cout << "the end of the word " << i + 1 <<" : " << ctr1 << endl;
        ctr1--;
        do
        {
            ctr1++;
        }
        while ( (str[ctr1] == ' ') || (str[ctr1] == '\t') || (str[ctr1] == '|') || (str[ctr1] == ',') || (str[ctr1] == '#') || (str[ctr1] == '<') || (str[ctr1] == '>') || (str[ctr1] == ')') || (str[ctr1] == '(') || (str[ctr1] == '"')  ) ;
//      cout << "the beginning of the word " << i + 2 <<" : " << ctr1 << endl;
        //skip all spaces
    }
    
    return ctr1;
}

int get_word_ending(int ctr1, char str[]) //looks for the end of the word starting with ctr1;
{
    do
    {
        ctr1++;
        if (str[ctr1] == 0)
        {
            return ctr1;
        }
    }
    while( (str[ctr1] != ' ') && (str[ctr1] != '\t') && (str[ctr1] != '|') && (str[ctr1] != ',') && (str[ctr1] != '#') && (str[ctr1] != '<') && (str[ctr1] != '>')  && (str[ctr1] != ')')  && (str[ctr1] != '(')  && (str[ctr1] != '"') );
    
    return ctr1;
}

int get_number(int ctr1, int ctr2, char str[]) // gets the word starting with ctr1 and ending with ctr2, translates it into a number if possible 
{
    char string[100];
    int i;
    for (i = ctr1; i < ctr2; i++)
    {
        string[i - ctr1] = str[i];
    }
    // cout << "i -ctr1 " << i - ctr1 << endl;
    string[i - ctr1] = 0;
    // cout << string << endl;
    return atoi(string);
}

void q_sort(std::vector<int>& numbers, std::vector<int>& addition, int left, int right)
{
    int pivot, pivot1, l_hold, r_hold;
    
    l_hold = left;
    r_hold = right;
    pivot = numbers[left];
    pivot1 = addition[left];
    while (left < right)
    {
        while ((numbers[right] >= pivot) && (left < right))
            right--;
        if (left != right)
        {
            numbers[left] = numbers[right];
            addition[left] = addition[right];
            left++;
        }
        while ((numbers[left] <= pivot) && (left < right))
            left++;
        if (left != right)
        {
            numbers[right] = numbers[left];
            addition[right] = addition[left];
            right--;
        }
    }
    numbers[left] = pivot;
    addition[left] = pivot1;
    pivot = left;
    left = l_hold;
    right = r_hold;
    if (left < pivot)
        q_sort(numbers, addition, left, pivot-1);
    if (right > pivot)
        q_sort(numbers, addition, pivot+1, right);
}

double incomplete_gamma(double a, double z)
{
    assert(a > -2);
    assert(z > 0);
    return std::exp(-z)*( -(std::exp(a*std::log(z))/a)   -  (std::exp((a+1)*std::log(z))/(a*(a+1)))) + (1./(a*(a+1))) * boost::math::tgamma(a+2, z);
}

double finite_avg_degree(int kmax, double kp_max, double avg_kp, double tpr, double gamma)
{
    return (kmax/kp_max) * (avg_kp / (tpr * finite_alpha(kp_max, gamma)));
}

double incomplete_negative_gamma_function(double x, double avg_k, double avg_kp, double tpr, double gamma, double k_obs, double kmax)
{
    double ff;
    ff = finite_avg_degree(kmax, x, avg_kp, tpr, gamma) * finite_alpha (x, gamma) * finite_alpha (x, gamma) * tpr;
    ff /=  (1 - (gamma - 1) * exp((gamma - 1)  * std::log(    (finite_avg_degree(kmax, x, avg_kp, tpr, gamma)/ avg_kp)  * tpr * finite_alpha(x, gamma)  ) ) * incomplete_gamma(1-gamma, (finite_avg_degree(kmax, x, avg_kp, tpr, gamma)/ avg_kp)  * tpr * finite_alpha(x, gamma)  ));
    return ff - k_obs;
}

void reverse(std::vector<int>& numbers)
{
    int i, size = numbers.size();
    std::vector<int> list(size);
    for (i = 0; i< size; i++)
    {
        list[i] = numbers[i];
    }
    for (i = 0; i< size; i++)
    {
        numbers[i] = list[size - i - 1];
    }
}