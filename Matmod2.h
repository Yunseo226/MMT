#include <bits/stdc++.h>
using namespace std;

class Matmod2{
public:
    int row;
    int column;
    vector<vector<bool>> mat;

    Matmod2(){
        ;
    }

    //make identity matrix
    Matmod2(int n){
        row = n;
        column = n;

        for(int i = 0; i < n; i++){
            vector<bool> vec(n);
            mat.push_back(vec);
        }

        for(int i = 0; i < n; i++){
            mat[i][i] = true;
        }
    }

    //make zero matrix
    Matmod2(int row, int column){
        this->row = row;
        this->column = column;

        for(int i = 0; i < row; i++){
            vector<bool> vec(column);
            mat.push_back(vec);
        }
    }

    //make column vector
    Matmod2(vector<int> y){
        row = y.size();
        column = 1;
        for(auto i : y){
            mat.push_back({(bool)i});
        }
    }

    //with vector<vector<bool>>
    Matmod2(vector<vector<bool>> _mat){
        row = _mat.size();
        column = _mat[0].size();
        mat = _mat;
    }

    //from text file to matrix.
    //only used for initializing function. The reason why this is not obvious function. (initially transpose)
    Matmod2(vector<string> data){
        this->column = data.size();
        this->row = data[0].size();
        for(int i = 0; i < row; i++){
            vector<bool> vec(column);
            mat.push_back(vec);
        }

        for(int i = 0; i < row; i++){
            for(int j = 0; j < column; j++){
                mat[i][j] = (data[j][i] == '1');
            }
        }        
    }

    Matmod2 operator*(Matmod2& ref){
        if(this->column != ref.row){
            cout << "error: matrix size not compatible for multiplication" << endl;
            return *this;
        }

        auto A = Matmod2(this->row, ref.column);
        for(int i = 0; i < A.row; i++){
            for(int j = 0; j < A.column; j++){
                for(int k = 0; k < this->column; k++){
                    A.mat[i][j] = (A.mat[i][j] != (this->mat[i][k]*ref.mat[k][j]));
                }
            }
        }

        return A;
    }

    Matmod2 operator+(Matmod2& ref){
         if(this->row != ref.row || this->column != ref.column){
            cout << "error: matrix size not compatible for addition" << endl;
            return *this;
        }

        auto A = Matmod2(this->row, this->column);
        for(int i = 0; i < A.row; i++){
            for(int j = 0; j < A.column; j++){
                A.mat[i][j] = (this->mat[i][j] != ref.mat[i][j]);
            }
        }

        return A;       
    }

    void transpose(){
        int tmp = row;
        row = column;
        column = tmp;

        vector<vector<bool>> new_mat;

        for(int i = 0; i < row; i++){
            vector<bool> vec(column);
            new_mat.push_back(vec);
        }

        for(int i = 0; i < row; i++){
            for(int j = 0; j < column; j++){
                new_mat[i][j] = mat[j][i];
            }
        }

        mat = new_mat;
    }

    void concat(Matmod2& ref){
        if(this->row != ref.row){
            cout << "error: matrix size not compatible for concatenation" << endl;
            return;
        }

        this->column += ref.column;
        for(int i = 0; i < this->row; i++){
            for(int j = 0; j < ref.column; j++){
                this->mat[i].push_back(ref.mat[i][j]);
            }
        }

    }

    void print_size(){
        cout << "row: " << row << ", column: " << column << endl;
    }

    void print(){
        for (int i=0; i<row; i++){
            for(int j = 0; j < column; j++){
                cout << mat[i][j] << " ";
            }
            cout << endl;
        }
    }

};

Matmod2 shuffle(Matmod2 A, vector<int> perm){
    auto B = A;
    B.transpose();
    auto C = B;
    for(int i = 0; i < perm.size(); i++){
        C.mat[perm[i]] = B.mat[i];
    }
    C.transpose();

    return C;
}

//cut row: a~b, column: x~y
Matmod2 cut(Matmod2 A, int a, int b, int x, int y){
    Matmod2 B(b-a+1, y-x+1);
    for(int i = 0; i < B.row; i++){
        for(int j = 0; j < B.column; j++){
            B.mat[i][j] = A.mat[a + i][x + j];
        }
    }

    return B;
}