#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include <cctype>

using namespace std;

string clean(string s){
    string ans;

    for (char c : s){
        if (c != '\r' && c != '\n'){
            ans.push_back(c);
        }
    }

    while (!ans.empty() && isspace(ans.front())){
        ans.erase(ans.begin());
    }

    while (!ans.empty() && isspace(ans.back())){
        ans.pop_back();
    }

    if (ans.size() >= 3 &&
        (unsigned char)ans[0] == 0xEF &&
        (unsigned char)ans[1] == 0xBB &&
        (unsigned char)ans[2] == 0xBF){
        ans = ans.substr(3);
    }

    return ans;
}

class DateUtil {
public:
    static double convert(string s){
        s = clean(s);

        char en = s.back();
        s.pop_back();

        double val = stod(s);

        if (en == 'D' || en == 'd') return val;
        if (en == 'W' || en == 'w') return val*7.0;
        if (en == 'M') return val*30.0;
        if (en == 'Y' || en == 'y') return val*360.0;

        throw runtime_error("Invalid duration: " + s + en);
    }
};

class MathUtil {
public:
    static vector<double> solve3(vector<vector<double>> a){
        int n = 3;

        for (int col = 0; col < n; col++){
            int piv = col;

            for (int row = col; row < n; row++){
                if (abs(a[row][col]) > abs(a[piv][col])){
                    piv = row;
                }
            }

            swap(a[col], a[piv]);

            double div = a[col][col];

            if (abs(div) < 1e-14){
                throw runtime_error("Singular matrix in solve3");
            }

            for (int j = col; j <= n; j++){
                a[col][j] /= div;
            }

            for (int row = 0; row < n; row++){
                if (row == col) continue;

                double val = a[row][col];

                for (int j = col; j <= n; j++){
                    a[row][j] -= val*a[col][j];
                }
            }
        }

        return {a[0][3], a[1][3], a[2][3]};
    }
};

class Interpolator {
public:
    virtual double get(vector<double> days, double duration, map<double,double> dfs) = 0;
    virtual ~Interpolator() {}
};

class LinearInterpolator : public Interpolator {
public:
    double get(vector<double> days, double duration, map<double,double> dfs) override {
        auto it = lower_bound(days.begin(), days.end(), duration);

        if (it != days.end() && abs(*it-duration) < 1e-9){
            return dfs[*it];
        }

        double t_prev, t_next;
        double df_prev, df_next;

        if (it == days.begin()){
            t_prev = 0.0;
            t_next = *it;
            df_prev = 1.0;
            df_next = dfs[t_next];
        }
        else if (it == days.end()){
            it--;
            return dfs[*it];
        }
        else{
            t_next = *it;
            df_next = dfs[t_next];

            it--;
            t_prev = *it;
            df_prev = dfs[t_prev];
        }

        double ln_df = log(df_prev) + 
            (log(df_next)-log(df_prev))*(duration-t_prev)/(t_next-t_prev);

        return exp(ln_df);
    }
};

class AveragedQuadraticInterpolator : public Interpolator {
private:
    LinearInterpolator lin;

    double quadratic_val(double x1, double x2, double x3,
                         double x, map<double,double> &dfs){
        double y1 = log(dfs[x1]);
        double y2 = log(dfs[x2]);
        double y3 = log(dfs[x3]);

        vector<vector<double>> a = {
            {x1*x1, x1, 1.0, y1},
            {x2*x2, x2, 1.0, y2},
            {x3*x3, x3, 1.0, y3}
        };

        vector<double> ans = MathUtil::solve3(a);

        return ans[0]*x*x + ans[1]*x + ans[2];
    }

public:
    double get(vector<double> days, double duration, map<double,double> dfs) override {
        int n = days.size();

        auto it = lower_bound(days.begin(), days.end(), duration);

        if (it != days.end() && abs(*it-duration) < 1e-9){
            return dfs[*it];
        }

        if (it == days.begin() || it == days.end()){
            return lin.get(days, duration, dfs);
        }

        int r = it - days.begin();
        int l = r - 1;

        if (l == 0){
            return lin.get(days, duration, dfs);
        }

        if (r == n-1){
            double val = quadratic_val(days[l-1], days[l], days[r], duration, dfs);
            return exp(val);
        }

        double val1 = quadratic_val(days[l-1], days[l], days[r], duration, dfs);
        double val2 = quadratic_val(days[l], days[r], days[r+1], duration, dfs);

        double w1 = (days[r] - duration)/(days[r] - days[l]);
        double w2 = (duration - days[l])/(days[r] - days[l]);

        return exp(w1*val1 + w2*val2);
    }
};

class DiscountCurve {
public:
    vector<double> days;
    map<double,double> dfs;

    void add_node(double t, double df){
        days.push_back(t);
        sort(days.begin(), days.end());
        dfs[t] = df;
    }

    double get_df(double t, Interpolator &interp){
        return interp.get(days, t, dfs);
    }
};

class CashCurveBuilder {
public:
    static double discount_factor(double days, double rate){
        return 1.0/(1.0 + rate*days/360.0);
    }

    DiscountCurve build(vector<double> days, vector<double> rates){
        DiscountCurve curve;

        for (int i = 0; i < (int)days.size(); i++){
            double df = discount_factor(days[i], rates[i]);
            curve.add_node(days[i], df);
        }

        return curve;
    }
};

class SwapCurveBuilder {
private:
    double swap_error(double x, double T, double S,
                      vector<double> days,
                      map<double,double> dfs,
                      Interpolator &interp){
        days.push_back(T);
        sort(days.begin(), days.end());
        dfs[T] = x;

        double fixed = 0.0;
        double prev = 0.0;

        for (double t = 180.0; t <= T + 1e-9; t += 180.0){
            double alpha = (t-prev)/360.0;
            double df = interp.get(days, t, dfs);

            fixed += S*alpha*df;
            prev = t;
        }

        double floating = 1.0 - x;

        return fixed - floating;
    }

    double solve_swap_df(double T, double S,
                         vector<double> days,
                         map<double,double> dfs,
                         Interpolator &interp){
        if (T <= 180.0){
            return CashCurveBuilder::discount_factor(T, S);
        }

        double lo = 1e-8;
        double hi = 2.0;

        for (int i = 0; i < 100; i++){
            double mid = (lo+hi)/2.0;
            double val = swap_error(mid, T, S, days, dfs, interp);

            if (val > 0){
                hi = mid;
            }
            else{
                lo = mid;
            }
        }

        return (lo+hi)/2.0;
    }

public:
    DiscountCurve build(vector<double> input_days,
                        vector<double> swap_rates,
                        Interpolator &interp){
        DiscountCurve curve;

        for (int i = 0; i < (int)input_days.size(); i++){
            double T = input_days[i];
            double S = swap_rates[i];

            double df = solve_swap_df(T, S, curve.days, curve.dfs, interp);
            curve.add_node(T, df);
        }

        return curve;
    }
};

class CSVReader {
public:
    static vector<vector<string>> read(string file){
        ifstream fin(file);

        if (!fin.is_open()){
            throw runtime_error("Could not open input file");
        }

        vector<vector<string>> rows;
        string line;

        while (getline(fin, line)){
            stringstream ss(line);
            string cell;
            vector<string> row;

            while (getline(ss, cell, ',')){
                row.push_back(clean(cell));
            }

            bool empty = true;
            for (auto &x : row){
                if (!x.empty()) empty = false;
            }

            if (!empty){
                rows.push_back(row);
            }
        }

        return rows;
    }
};

struct SwapResult {
    double pv;
    double par_rate;
};

class SwapPricer {
public:
    static double freq_to_days(string s){
        s = clean(s);

        char en = s.back();
        s.pop_back();

        double val = stod(s);

        if (en == 'm' || en == 'M') return val*30.0;
        if (en == 'y' || en == 'Y') return val*360.0;

        throw runtime_error("Invalid frequency: " + s + en);
    }

    static SwapResult price(double notional,
                            double fixed_rate,
                            double maturity,
                            double fixed_freq,
                            DiscountCurve &curve,
                            Interpolator &interp){
        double annuity = 0.0;
        double prev = 0.0;

        for (double t = fixed_freq; t <= maturity + 1e-9; t += fixed_freq){
            double alpha = (t-prev)/360.0;
            double df = curve.get_df(t, interp);

            annuity += alpha*df;
            prev = t;
        }

        double dfT = curve.get_df(maturity, interp);

        double pv_fixed = notional*fixed_rate*annuity;
        double pv_float = notional*(1.0 - dfT);

        SwapResult res;
        res.pv = pv_float - pv_fixed;
        res.par_rate = (1.0 - dfT)/annuity;

        return res;
    }
};


class RiskCalculator {
public:
    static double dnode_dcash(double t, double c){
        double val = t/360.0;
        return -val/((1.0 + c*val)*(1.0 + c*val));
    }

    static double lagrange(double x, double a, double b, double c, int pos){
        if (pos == 0){
            return ((x-b)*(x-c))/((a-b)*(a-c));
        }
        if (pos == 1){
            return ((x-a)*(x-c))/((b-a)*(b-c));
        }

        return ((x-a)*(x-b))/((c-a)*(c-b));
    }

    static double dquad_dnode(double t,
                              int i,
                              double a,
                              double b,
                              double c,
                              DiscountCurve &curve,
                              AveragedQuadraticInterpolator &aq){
        double mi = curve.days[i];

        if (abs(mi-a) > 1e-9 && abs(mi-b) > 1e-9 && abs(mi-c) > 1e-9){
            return 0.0;
        }

        int pos = 0;

        if (abs(mi-b) < 1e-9){
            pos = 1;
        }
        else if (abs(mi-c) < 1e-9){
            pos = 2;
        }

        double L = lagrange(t, a, b, c, pos);
        double df = curve.get_df(t, aq);
        double dfi = curve.dfs[mi];

        return df*L/dfi;
    }

    static double dlinear_dnode(double t,
                               int i,
                               vector<double> days,
                               DiscountCurve &curve,
                               LinearInterpolator &lin){
        int n = days.size();
        double mi = days[i];
        double dfi = curve.dfs[mi];
        double df = curve.get_df(t, lin);

        if (abs(t-mi) < 1e-9){
            return 1.0;
        }

        if (i == 0){
            if (0.0 <= t && t <= days[0]){
                double w = t/days[0];
                return w*df/dfi;
            }
        }
        else{
            if (days[i-1] <= t && t <= days[i]){
                double w = (t-days[i-1])/(days[i]-days[i-1]);
                return w*df/dfi;
            }
        }

        if (i+1 < n && days[i] <= t && t <= days[i+1]){
            double w = (days[i+1]-t)/(days[i+1]-days[i]);
            return w*df/dfi;
        }

        return 0.0;
    }

    static double daq_dnode(double t,
                            int i,
                            vector<double> days,
                            DiscountCurve &curve,
                            AveragedQuadraticInterpolator &aq,
                            LinearInterpolator &lin){
        int n = days.size();

        if (t <= days[1] || t > days[n-1]){
            return dlinear_dnode(t, i, days, curve, lin);
        }

        auto it = lower_bound(days.begin(), days.end(), t);

        if (it != days.end() && abs(*it-t) < 1e-9){
            if (abs(days[i]-t) < 1e-9){
                return 1.0;
            }
            return 0.0;
        }

        if (it == days.begin() || it == days.end()){
            return dlinear_dnode(t, i, days, curve, lin);
        }

        int r = it - days.begin();
        int l = r - 1;

        if (l == 0){
            return dlinear_dnode(t, i, days, curve, lin);
        }

        if (r == n-1){
            return dquad_dnode(t, i, days[l-1], days[l], days[r], curve, aq);
        }

        double w1 = (days[r] - t)/(days[r] - days[l]);
        double w2 = (t - days[l])/(days[r] - days[l]);

        double val1 = dquad_dnode(t, i, days[l-1], days[l], days[r], curve, aq);
        double val2 = dquad_dnode(t, i, days[l], days[r], days[r+1], curve, aq);

        return w1*val1 + w2*val2;
    }

    static vector<double> cash_linear(double notional,
                                      double fixed_rate,
                                      double maturity,
                                      double fixed_freq,
                                      vector<double> days,
                                      vector<double> cash_rates,
                                      DiscountCurve &curve,
                                      LinearInterpolator &lin){
        int n = days.size();
        vector<double> ans(n, 0.0);

        for (int i = 0; i < n; i++){
            double val = 0.0;
            double prev = 0.0;

            val += -notional*dlinear_dnode(maturity, i, days, curve, lin);

            for (double t = fixed_freq; t <= maturity + 1e-9; t += fixed_freq){
                double alpha = (t-prev)/360.0;
                val += -notional*fixed_rate*alpha*dlinear_dnode(t, i, days, curve, lin);
                prev = t;
            }

            ans[i] = val*dnode_dcash(days[i], cash_rates[i]);
        }

        return ans;
    }

    static vector<double> cash_aq(double notional,
                                  double fixed_rate,
                                  double maturity,
                                  double fixed_freq,
                                  vector<double> days,
                                  vector<double> cash_rates,
                                  DiscountCurve &curve,
                                  AveragedQuadraticInterpolator &aq,
                                  LinearInterpolator &lin){
        int n = days.size();
        vector<double> ans(n, 0.0);

        for (int i = 0; i < n; i++){
            double val = 0.0;
            double prev = 0.0;

            val += -notional*daq_dnode(maturity, i, days, curve, aq, lin);

            for (double t = fixed_freq; t <= maturity + 1e-9; t += fixed_freq){
                double alpha = (t-prev)/360.0;
                val += -notional*fixed_rate*alpha*daq_dnode(t, i, days, curve, aq, lin);
                prev = t;
            }

            ans[i] = val*dnode_dcash(days[i], cash_rates[i]);
        }

        return ans;
    }

    static double swap_eq_derivative_node_linear(double T,
                                                 double S,
                                                 int node,
                                                 vector<double> days,
                                                 DiscountCurve &curve,
                                                 LinearInterpolator &lin){
        double val = dlinear_dnode(T, node, days, curve, lin);
        double prev = 0.0;

        for (double t = 180.0; t <= T + 1e-9; t += 180.0){
            double alpha = (t-prev)/360.0;
            val += S*alpha*dlinear_dnode(t, node, days, curve, lin);
            prev = t;
        }

        return val;
    }

    static double swap_eq_derivative_node_aq(double T,
                                             double S,
                                             int node,
                                             vector<double> days,
                                             DiscountCurve &curve,
                                             AveragedQuadraticInterpolator &aq,
                                             LinearInterpolator &lin){
        double val = daq_dnode(T, node, days, curve, aq, lin);
        double prev = 0.0;

        for (double t = 180.0; t <= T + 1e-9; t += 180.0){
            double alpha = (t-prev)/360.0;
            val += S*alpha*daq_dnode(t, node, days, curve, aq, lin);
            prev = t;
        }

        return val;
    }

    static double calib_annuity_linear(double T,
                                       double S,
                                       vector<double> days,
                                       DiscountCurve &curve,
                                       LinearInterpolator &lin){
        double ans = 0.0;
        double prev = 0.0;

        for (double t = 180.0; t <= T + 1e-9; t += 180.0){
            double alpha = (t-prev)/360.0;
            ans += alpha*curve.get_df(t, lin);
            prev = t;
        }

        return ans;
    }

    static double calib_annuity_aq(double T,
                                   double S,
                                   vector<double> days,
                                   DiscountCurve &curve,
                                   AveragedQuadraticInterpolator &aq){
        double ans = 0.0;
        double prev = 0.0;

        for (double t = 180.0; t <= T + 1e-9; t += 180.0){
            double alpha = (t-prev)/360.0;
            ans += alpha*curve.get_df(t, aq);
            prev = t;
        }

        return ans;
    }

    static vector<vector<double>> swap_node_sens_linear(vector<double> days,
                                                        vector<double> swap_rates,
                                                        DiscountCurve &curve,
                                                        LinearInterpolator &lin){
        int n = days.size();
        vector<vector<double>> sens(n, vector<double>(n, 0.0));

        for (int i = 0; i < n; i++){
            double T = days[i];
            double S = swap_rates[i];

            if (T <= 180.0){
                double val = T/360.0;
                sens[i][i] = -val/((1.0 + S*val)*(1.0 + S*val));
                continue;
            }

            DiscountCurve cur;
            for (int p = 0; p <= i; p++){
                cur.add_node(days[p], curve.dfs[days[p]]);
            }

            vector<double> cur_days(days.begin(), days.begin()+i+1);

            double ann = calib_annuity_linear(T, S, cur_days, cur, lin);
            double den = swap_eq_derivative_node_linear(T, S, i, cur_days, cur, lin);

            for (int k = 0; k <= i; k++){
                double val = (k == i ? ann : 0.0);

                for (int p = 0; p < i; p++){
                    double der = swap_eq_derivative_node_linear(T, S, p, cur_days, cur, lin);
                    val += der*sens[p][k];
                }

                sens[i][k] = -val/den;
            }
        }

        return sens;
    }

    static vector<vector<double>> swap_node_sens_aq(vector<double> days,
                                                    vector<double> swap_rates,
                                                    DiscountCurve &curve,
                                                    AveragedQuadraticInterpolator &aq,
                                                    LinearInterpolator &lin){
        int n = days.size();
        vector<vector<double>> sens(n, vector<double>(n, 0.0));

        for (int i = 0; i < n; i++){
            double T = days[i];
            double S = swap_rates[i];

            if (T <= 180.0){
                double val = T/360.0;
                sens[i][i] = -val/((1.0 + S*val)*(1.0 + S*val));
                continue;
            }

            DiscountCurve cur;
            for (int p = 0; p <= i; p++){
                cur.add_node(days[p], curve.dfs[days[p]]);
            }

            vector<double> cur_days(days.begin(), days.begin()+i+1);

            double ann = calib_annuity_aq(T, S, cur_days, cur, aq);
            double den = swap_eq_derivative_node_aq(T, S, i, cur_days, cur, aq, lin);

            for (int k = 0; k <= i; k++){
                double val = (k == i ? ann : 0.0);

                for (int p = 0; p < i; p++){
                    double der = swap_eq_derivative_node_aq(T, S, p, cur_days, cur, aq, lin);
                    val += der*sens[p][k];
                }

                sens[i][k] = -val/den;
            }
        }

        return sens;
    }

    static vector<double> swap_linear(double notional,
                                      double fixed_rate,
                                      double maturity,
                                      double fixed_freq,
                                      vector<double> days,
                                      vector<double> swap_rates,
                                      DiscountCurve &curve,
                                      LinearInterpolator &lin){
        int n = days.size();
        vector<double> node_val(n, 0.0);
        vector<double> ans(n, 0.0);

        for (int i = 0; i < n; i++){
            double val = 0.0;
            double prev = 0.0;

            val += -notional*dlinear_dnode(maturity, i, days, curve, lin);

            for (double t = fixed_freq; t <= maturity + 1e-9; t += fixed_freq){
                double alpha = (t-prev)/360.0;
                val += -notional*fixed_rate*alpha*dlinear_dnode(t, i, days, curve, lin);
                prev = t;
            }

            node_val[i] = val;
        }

        vector<vector<double>> sens = swap_node_sens_linear(days, swap_rates, curve, lin);

        for (int k = 0; k < n; k++){
            for (int i = k; i < n; i++){
                ans[k] += node_val[i]*sens[i][k];
            }
        }

        return ans;
    }

    static vector<double> swap_aq(double notional,
                                  double fixed_rate,
                                  double maturity,
                                  double fixed_freq,
                                  vector<double> days,
                                  vector<double> swap_rates,
                                  DiscountCurve &curve,
                                  AveragedQuadraticInterpolator &aq,
                                  LinearInterpolator &lin){
        int n = days.size();
        vector<double> node_val(n, 0.0);
        vector<double> ans(n, 0.0);

        for (int i = 0; i < n; i++){
            double val = 0.0;
            double prev = 0.0;

            val += -notional*daq_dnode(maturity, i, days, curve, aq, lin);

            for (double t = fixed_freq; t <= maturity + 1e-9; t += fixed_freq){
                double alpha = (t-prev)/360.0;
                val += -notional*fixed_rate*alpha*daq_dnode(t, i, days, curve, aq, lin);
                prev = t;
            }

            node_val[i] = val;
        }

        vector<vector<double>> sens = swap_node_sens_aq(days, swap_rates, curve, aq, lin);

        for (int k = 0; k < n; k++){
            for (int i = k; i < n; i++){
                ans[k] += node_val[i]*sens[i][k];
            }
        }

        return ans;
    }
};

class CSVWriter {
public:
    static void write_output(string file,
                             double q1a, double q1b, double q1c, double q1d,
                             SwapResult q21a, SwapResult q21b,
                             SwapResult q21c, SwapResult q21d,
                             vector<double> q22a,
                             vector<double> q22b,
                             vector<double> q22c,
                             vector<double> q22d){
        ofstream fout(file);

        if (!fout.is_open()){
            throw runtime_error("Could not open output file");
        }

        fout << fixed << setprecision(10);

        fout << q1a << "," << q1b << "," << q1c << "," << q1d << "\n";

        fout << q21a.pv << "," << q21b.pv << ","
             << q21c.pv << "," << q21d.pv << "\n";

        fout << q21a.par_rate << "," << q21b.par_rate << ","
             << q21c.par_rate << "," << q21d.par_rate << "\n";

        for (int i = 0; i < (int)q22a.size(); i++){
            fout << q22a[i] << "," << q22b[i] << "," << q22c[i] << "," << q22d[i] << "\n";
        }
    }
};
bool is_number(string s){
    s = clean(s);

    if (s.empty()) return false;

    int i = 0;

    if (s[0] == '-' || s[0] == '+'){
        i++;
    }

    bool seen = false;

    for (; i < (int)s.size(); i++){
        if (!isdigit(s[i])){
            return false;
        }
        seen = true;
    }

    return seen;
}

int main(){
    try {
        vector<vector<string>> rows = CSVReader::read("Input.csv");

        int start = 0;

        while (start < (int)rows.size() && 
               (rows[start].empty() || !is_number(rows[start][0]))){
            start++;
        }

        if (start == (int)rows.size()){
            throw runtime_error("Could not find maturity count row");
        }

        int n = stoi(rows[start][0]);

        vector<double> days;
        vector<double> cash_rates;
        vector<double> swap_rates;

        for (int i = start+1; i <= start+n; i++){
            if ((int)rows[i].size() < 3){
                throw runtime_error("Invalid market data row");
            }

            double t = DateUtil::convert(rows[i][0]);
            double cash = stod(rows[i][1])/100.0;
            double swap = stod(rows[i][2])/100.0;

            days.push_back(t);
            cash_rates.push_back(cash);
            swap_rates.push_back(swap);
        }

        double query_time = stod(rows[start+n+1][0]);

        LinearInterpolator lin;
        AveragedQuadraticInterpolator aq;

        CashCurveBuilder cash_builder;
        SwapCurveBuilder swap_builder;

        DiscountCurve cash_curve = cash_builder.build(days, cash_rates);

        DiscountCurve swap_curve_lin = swap_builder.build(days, swap_rates, lin);
        DiscountCurve swap_curve_aq = swap_builder.build(days, swap_rates, aq);

        double q1a = cash_curve.get_df(query_time, lin);
        double q1b = cash_curve.get_df(query_time, aq);
        double q1c = swap_curve_lin.get_df(query_time, lin);
        double q1d = swap_curve_aq.get_df(query_time, aq);

        int swap_row = start+n+2;

        double notional = 100.0;
        double fixed_rate = stod(rows[swap_row][0])/100.0;
        double swap_maturity = DateUtil::convert(rows[swap_row][1]);
        double fixed_freq = SwapPricer::freq_to_days(rows[swap_row][2]);

        SwapResult q21a = SwapPricer::price(notional, fixed_rate, swap_maturity, fixed_freq, cash_curve, lin);
        SwapResult q21b = SwapPricer::price(notional, fixed_rate, swap_maturity, fixed_freq, cash_curve, aq);
        SwapResult q21c = SwapPricer::price(notional, fixed_rate, swap_maturity, fixed_freq, swap_curve_lin, lin);
        SwapResult q21d = SwapPricer::price(notional, fixed_rate, swap_maturity, fixed_freq, swap_curve_aq, aq);

        vector<double> q22a = RiskCalculator::cash_linear(notional,
                                                           fixed_rate,
                                                           swap_maturity,
                                                           fixed_freq,
                                                           days,
                                                           cash_rates,
                                                           cash_curve,
                                                           lin);

        vector<double> q22b = RiskCalculator::cash_aq(notional,
                                                     fixed_rate,
                                                     swap_maturity,
                                                     fixed_freq,
                                                     days,
                                                     cash_rates,
                                                     cash_curve,
                                                     aq,
                                                     lin);

        vector<double> q22c = RiskCalculator::swap_linear(notional,
                                                          fixed_rate,
                                                          swap_maturity,
                                                          fixed_freq,
                                                          days,
                                                          swap_rates,
                                                          swap_curve_lin,
                                                          lin);

        vector<double> q22d = RiskCalculator::swap_aq(notional,
                                                      fixed_rate,
                                                      swap_maturity,
                                                      fixed_freq,
                                                      days,
                                                      swap_rates,
                                                      swap_curve_aq,
                                                      aq,
                                                      lin);

        CSVWriter::write_output("Output.csv",
                                q1a, q1b, q1c, q1d,
                                q21a, q21b, q21c, q21d,
                                q22a, q22b, q22c, q22d);

    }
    catch (const exception &e){
        cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}