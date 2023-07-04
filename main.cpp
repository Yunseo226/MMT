#include "Matmod2.h"

int n;
int k;
int w;
Matmod2 H;
Matmod2 s;

bool init(){
    string filename("input.txt");
    vector<string> lines;
    string line;

    ifstream input_file(filename);
    if (!input_file.is_open()) {
        cerr << "Could not open the file - '"
             << filename << "'" << endl;
        return false;
    }

    while (getline(input_file, line)){
        lines.push_back(line);
    }

    input_file.close();

    n = stoi(lines[1]); 
    w = stoi(lines[5]); 
    k = n - lines[7].size();

    char c;
    int i = 7;
    vector<string> H_str;
    vector<string> s_str;
    while(true){
        H_str.push_back(lines[i]);
        i++;
        c = lines[i][0];
        if(c == '#') break;
    }

    s_str.push_back(lines[i+1]);
    s = Matmod2(s_str);

    auto extra = Matmod2(H_str);
    H = Matmod2(n-k);
    H.concat(extra);

    return true;    
}

//gauss with only first p columns, and check their rank.
bool gauss(Matmod2& A, Matmod2& b, int p){
    if(p > n-k) {
        cout << "error: two big integer" << endl;
        return false;
    }

    bool is_zero = true;
    for(int t = 0; t < p; t++){
        for(int i = 0; i < n - k; i++){
            if(A.mat[i][t]){
                swap(A.mat[i], A.mat[t]);
                swap(b.mat[i], b.mat[t]);
                for(int j = i + 1; j < n - k; j++){
                    if(A.mat[j][t]){
                        for(int x = 0; x < n; x++){
                            A.mat[j][x] = (A.mat[j][x] != A.mat[t][x]);
                        }
                        b.mat[j][0] = (b.mat[j][0] != b.mat[t][0]);
                    }
                }
                is_zero = false;
                break;
            }
        }
        if(is_zero){
            return false;
        }
    }

    return true;
}
 
int main()
{
    /* input matrix */
    if(!init()){
        cout << "error: cannot initialize" << endl;
        return 0;
    }
 
    int p = 4;   //choose between 0 <= p <= W, (should it be multiple of 4)
    int l = 6;   //choose between 0 <= l <= n - k, should be k == l (mod 2)

    int d[31][31];
    for (int i = 0; i < 31; i++) {
        d[i][0] = 1;
    }
    for (int i = 1; i < 31; i++) {
        for (int j = 1; j <= i; j++) {
            d[i][j] = d[i - 1][j - 1] + d[i - 1][j];
        }
    } //pascal's triangle for combination

    int p1 = p/2;
    int p2 = p1/2;
    int l1 = (int)log2(d[p-1][p1-1]);

    vector<int> perm;
    for(int i = 0; i < n; i++){
        perm.push_back(i);
    }

    do {
        auto H_bar = shuffle(H, perm);

        if(!gauss(H_bar, s, n - k - l)){
            continue;
        }

        auto H1 = cut(H_bar, 0, n-k-l-1, n-k-l, n-1);
        auto H2 = cut(H_bar, n-k-l, n-1, n-k-l, n-1);
        auto s1 = cut(s, 0, n-k-l-1, 0, 0);
        auto s2 = cut(s, n-k-l, n-1, 0, 0);

        vector<pair<Matmod2, Matmod2>> L1, L2, L3, L4, LL1, LL2, L;
        vector<Matmod2> y1_set, y2_set;

        vector<int> ones((k+l)/2 - p2);
        for(int i = 0; i < p2; i++){
            ones.push_back(1);
        }
        vector<int> zero((k+l)/2);

        do{
            vector<int> y1 = ones;
            vector<int> y2 = zero;
            y1.insert(y1.end(), zero.begin(), zero.end());
            y2.insert(y2.end(), ones.begin(), ones.end());
            y1_set.push_back(Matmod2(y1));
            y2_set.push_back(Matmod2(y2));
        }
        while(next_permutation(y.begin(), y.end()));

        for(auto y1 : y1_set){
            L1.push_back(make_pair(y1, H1*y1));
            L3.push_back(make_pair(y1, H1*y1));
        }

        for(auto y2 : y2_set){
            L2.push_back(make_pair(y2, H2*y2));
            L4.push_back(make_pair(y2, H2*y2 + s2));
        }

    } while(next_permutation(perm.begin(), perm.end()));
 
    return 0;
}