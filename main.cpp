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
    H.concat_in_row(extra);

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
 
    int l;   //choose between 0 <= l <= n - k, (it must be k == l (mod 2))
    int p;   //choose between 0 <= p <=  min(w, k + l), (it must be multiple of 4)

    cout << "choose l between 0 and " << n-k << endl;
    cin >> l;
    cout << "choose p between 0 and " << (w > (k + l) ? k+l : w) << endl;
    cin >> p;

    p >> 2;
    p << 2;
    if(p == 0) p = 4;
    if(w == p) p-=4;

    l += (k-l)%2;
    if(l >= n-k) l-=2;
    if(l < 0) l+=2;

    auto start = chrono::steady_clock::now();


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
    int l1 = (int)log2(d[p][p1]);
    Matmod2 ans;

    vector<int> perm;
    for(int i = 0; i < n; i++){
        perm.push_back(i);
    }

    int iter = 0;

    do {
        auto H_bar = shuffle_col(H, perm);

        if(!gauss(H_bar, s, n - k - l)){
            continue;
        }

        auto H1 = cut(H_bar, 0, n-k-l-1, n-k-l, n-1);
        auto H2 = cut(H_bar, n-k-l, n-k-1, n-k-l, n-1);
        auto s1 = cut(s, 0, n-k-l-1, 0, 0);
        auto s2 = cut(s, n-k-l, n-k-1, 0, 0);

        vector<Matmod2> L;
        vector<pair<Matmod2, Matmod2>> L1, L2, L3, L4, LL1, LL2;
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
        while(next_permutation(ones.begin(), ones.end()));

        for(auto y1 : y1_set){
            L1.push_back(make_pair(y1, H2*y1));
            L3.push_back(make_pair(y1, H2*y1));
        }

        for(auto y2 : y2_set){
            L2.push_back(make_pair(y2, H2*y2));
            L4.push_back(make_pair(y2, H2*y2 + s2));
        }

        auto t = gen_randvec(l1);

        for(auto v: L1){
            for(auto w : L2){
                if(proj(v.second, l1) == proj(w.second, l1) + t){
                    LL1.push_back(make_pair(v.first + w.first, v.second + w.second));
                }
            }
        }

        for(auto v: L3){
            for(auto w : L4){
                if(proj(v.second, l1) == proj(w.second, l1) + t){
                    LL2.push_back(make_pair(v.first + w.first, v.second + w.second));
                }
            }
        }

        for(auto x1 : LL1){
            for(auto x2 : LL2){
                if(x1.second == x2.second){
                    L.push_back(x1.first + x2.first);
                }
            }
        }

        bool solved = false;
        for(auto e2 : L){
            auto e1 = H1*e2 + s1;
            if(wt(e1) <= w-p){
                e1.concat_in_col(e2);
                ans = shuffle_row(e1, perm);
                solved = true;
                break;
            }
        }

        if(solved) break;
        iter++;

    } while(next_permutation(perm.begin(), perm.end()));
 

    cout << "answer: ";
    for(int i = 0; i < ans.row; i++){
        cout << ans.mat[i][0];
    }
    cout << endl;

    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << "number of iteration: " << iter << endl;
    cout << "execution time: " << chrono::duration <double, milli> (diff).count() << " ms" << endl;

    return 0;
}