#include <iostream>
#include <vector>
#include <math.h> 
using namespace std;

class Bond{
    public:

    Bond(double F, double issue_date, int num_periods, int freq, vector<double> &c){//Constructor; hence, the similtude of class name.
        if(F >= 0) {_Face = F;};
        if(freq >= 1) { _cpnFreq = freq;};
        if(num_periods >= 1) { _numCpnPeriods= num_periods;};
        _issue = issue_date;

        _maturity = _issue + (double)_numCpnPeriods/(double)_cpnFreq; 
        //date of issue  + how many coupons you will recieve divided how many per year

       
        _cpnDate.resize(_numCpnPeriods);
        _cpnAmt.resize(_numCpnPeriods);

        setCoupons(c);

        for(int i= _issue; i<_numCpnPeriods-1; i++ ){
            test = 0.0;
            test = _issue+(_numCpnPeriods/_cpnFreq);
            _cpnDate[i] = test;
            cout<<test<<"values\n";
            if(i == _numCpnPeriods-1 && _maturity==_cpnDate[i]) cout<<"Probably correct\n"; 
        }
    };
    ~Bond( ){
            
    };
    

    void setFlatCoupons(double c){
        if(c < 0.0) c = 0.0;
        fill(_cpnAmt.begin(), _cpnAmt.end(), c);
    };
    void setCoupons(vector<double> &c){
        double temp;
        int size = c.size();
        int tempC = _numCpnPeriods;
        for(int i= 0; i< _numCpnPeriods; i++){
            if(!(c[i] >= 0)){ c[i]=0;};
            _cpnAmt.insert (_cpnAmt.begin()+i, c[i] );
            if( tempC> size ){ 
                temp = c.back();
                 _cpnAmt[size++] = temp;
                 tempC--; }
        }
    };
    double FairValue(double t0, double y) const
    {
        double B = 0;
        double dummy1 = 0;
        double dummy2 = 0;
        FV_duration(t0, y, B, dummy1, dummy2);
        return B;
    };
    double maturity() const { return _maturity; };
    double issue() const { return _issue; };

    int FV_duration(double t0, double y, double &B, double &Mac_dur, double &mod_dur) const 
    {
        B = 0, Mac_dur = 0, mod_dur = 0;
        double yD = (0.01 * y);
        const double tol = 1.0e-6;

        if(t0< _issue || t0 >= _maturity) return 1;

        for(int i = 0; i< _numCpnPeriods; i++){
            if(_cpnDate[i] >= t0 + tol ){
            B= B + (i/_cpnFreq)/(1+(yD/_cpnFreq))*(exp(_cpnDate[i]-t0));
            if(i == _numCpnPeriods-1){
               B= B+ _Face+(i/_cpnFreq)/(1+(yD/_cpnFreq))*(exp(_cpnDate[i]-t0));
            }

            
            
        }
        Mac_dur= Mac_dur+((1/B)*(_cpnDate[i]-t0)*(_cpnDate[i]/_cpnFreq)*(1+yD/_cpnFreq));
        }


        return 0;
    };
      private:
    double _Face;
    double _issue;
    double _maturity;
    int _cpnFreq;
    int _numCpnPeriods;
    vector<double> _cpnAmt;
    vector<double> _cpnDate; 
    double test;
};


int yield(double &y, int &num_iter , const Bond &bond, double B_target, double t0,
double tol=1.0e-4, int max_iter=100){
 num_iter = 0;
 y=y +0.0; 
if(B_target<0.0 || t0 < bond.issue() || t0 >= bond.maturity()) return 1;

double y_low = 0.0;
double y_high = 100.0;
double B_y_low = bond.FairValue(t0, y_low);
double diff_B_y_low = B_y_low - B_target;

if(abs(diff_B_y_low) <= tol) y=y_low; return 0;
y_high = 100.0; 
double B_y_high = bond.FairValue(t0, y_high);

double diff_B_y_high = B_y_high - B_target;

if(abs(diff_B_y_high) <= tol) y=y_high; return 0;



bool oppositeSigns =((diff_B_y_low * diff_B_y_high) > 0.0);
    if(oppositeSigns == true){
        y=0;
        return 1;
    }

    for (num_iter = 1; num_iter < max_iter; ++num_iter){
        y=(y_low + y_high)/2;
       double  B = bond.FairValue(t0, y);
        double diff_B = B -B_target;
        if(abs(diff_B) <= tol) return 0;

bool oppositeSigns =((diff_B_y_low * diff_B) > 0.0);
        if(oppositeSigns == true) y_low = y;
        else(y_high = y);

        if(abs(y_high - y_low) <= tol) return 0;



    }
    if(num_iter >= max_iter){
        y =0;
        return 1;
    }
};

// int main(){

// vector<double> _cpn;
// int Face;
// int freq;
// int nP;
// int i =0;
// string input;
// cout << "Enter a face value\n" <<endl;
// cin >> Face;
// cout << "Enter number of periods\n" <<endl;
// cin >> nP;
// cout <<"Enter frequency of coupon\n" <<endl;
// cin >>  freq;
// cout << "Enter coupon amounts and if all is entered, input stop\n" <<endl;
// while(cin >> input ){
//     double value = 0.0;
//     if(input == "stop") return 1;
    
//      value = stod(input);
//     cout << value<<"\n"<<endl;
//     i++;
    
    
//     _cpn.insert(_cpn.begin()+i, value);
// }
// Bond James(100, 0, nP=0, freq= 0, _cpn );




//     return 0;
// }


//Maturity() and issue() is a const to avoid modification of a date "set in stone"
int main()
{
  double F = 100.0;
  int f = 2;
  std::vector<double> c(1, 4.0);
  int num_cpn1 = 4;
  int num_cpn2 = 8;
  double issue1 = 0.0;
  double issue2 = -0.56;
  double t0 = 0.0;
  double y = 5.0;
  double FV = 0.0;
  double Mac_dur = 0.0;
  double mod_dur = 0.0;
   
  Bond bond1(F, issue1, num_cpn1, f, c);
  Bond bond2(F, issue2, num_cpn2, f, c);
  
  bond1.FV_duration(t0, y, FV, Mac_dur, mod_dur);
  std::cout << "#1a: " << FV << "   " << Mac_dur << "   " << mod_dur << std::endl;
    
  bond2.FV_duration(t0, y, FV, Mac_dur, mod_dur);
  std::cout << "#2a: " << FV << "   " << Mac_dur << "   " << mod_dur << std::endl;
  
  std::vector<double> cpn(100, 4.0);
  bond1.setCoupons(cpn);
  bond2.setCoupons(cpn);
  
  bond1.FV_duration(t0, y, FV, Mac_dur, mod_dur);
  std::cout << "#1b: " << FV << "   " << Mac_dur << "   " << mod_dur << std::endl;
  
  bond2.FV_duration(t0, y, FV, Mac_dur, mod_dur);
  std::cout << "#2b: " << FV << "   " << Mac_dur << "   " << mod_dur << std::endl;
  
  return 0;
}

// int main()
// {
//      double F = 100.0;
//   int f = 2;
//   std::vector<double> c(1, 4.0);
//   int num_cpn1 = 8;
//   double issue1 = 0.0;
//   double t0 = 0.73410543;
//   double y = 5.0;
  
//      Bond bond1(F, issue1, num_cpn1, f, c);
//     double yV= 109.08;
//     int iV = 0;
//    int y0= yield(yV, iV, bond1, 100.23410543 , 0.73410543 );
//     cout<<y0<<endl;
// }