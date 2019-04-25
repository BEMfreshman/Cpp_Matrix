//
// Created by YangWenlu on 2019/4/5.
//

#ifndef CPP_MATRIX_UTILITY_H
#define CPP_MATRIX_UTILITY_H

template <typename T>
class Utility
{
public:
    static const Matrix<T> ToRowEchelonForm(const Matrix<T>& mat, std::vector<size_t>& RowReturn,
                                      std::vector<size_t>& ColReturn, std::vector<double>& NumReturn, size_t& FirstTransFormTimes)
    {
        //采用消元将一个矩阵变为行最简矩阵
        //RowReturn，ColReturn，NumReturn返回的是消元时对哪一行哪一列加上了哪个数字将这一列的非主元项目削减为0（在LU分解中用于L矩阵的产生）
        //FirstTransFormTimes表示做了几次第一类变换（行变换）

        FirstTransFormTimes = 0; //初始化为0
        size_t Row = mat.GetNumRow();
        size_t Col = mat.GetNumCol();
        Matrix<T> matOutput;
        //输出行最简矩阵
        matOutput = mat;

        size_t row_counter = 0;
        size_t col_counter = 0;
        T PivotElement;

        for (col_counter = 0; col_counter < Col; col_counter++)
        {
            row_counter = col_counter;


            //对角线主元是否为0？如果是则交换，如果不是则继续
            //如果需要交换，则采用初等变换进行

            PivotElement = matOutput(row_counter, row_counter);
            size_t tmp = 0;
            while (abs(PivotElement) < EPS)
            {
                // 如果主元为0则在该列一直寻找，直到找到一个主元不为0的数字
                tmp++;

                if (row_counter + tmp == Row)
                {
                    //表示本列中所有元素均为0
                    break;
                }
                PivotElement = matOutput(row_counter + tmp, row_counter);
            }

            if (tmp != 0 && row_counter + tmp < Row)
            {
                //该行的对角主元确实为0
                //且在该列中确实发现不为0的元素

                //进行交换
                matOutput.FirstTypeTransForm(row_counter + tmp, row_counter);
                FirstTransFormTimes++;
            }

            if (row_counter + tmp == Row)
            {
                continue;
            }

            PivotElement = matOutput(row_counter, row_counter);



            for (size_t tmp_row_counter = row_counter + 1; tmp_row_counter < Col; tmp_row_counter++)
            {
                //从对角的下一行开始
                T num = matOutput(tmp_row_counter, col_counter);
                double times = num / PivotElement;
                // 倍数

                matOutput.ThirdTypeTransForm(row_counter, tmp_row_counter, -times);
                RowReturn.push_back(tmp_row_counter);
                ColReturn.push_back(col_counter);
                NumReturn.push_back(-times);
            }
        }
        return matOutput;
    }


    static Matrix<T>& ToRowEchelonForm(Matrix<T>& mat, std::vector<size_t>& RowReturn,
                                       std::vector<size_t>& ColReturn, std::vector<double>& NumReturn, size_t& FirstTransFormTimes)
    {
        //采用消元将一个矩阵变为行最简矩阵
        //RowReturn，ColReturn，NumReturn返回的是消元时对哪一行哪一列加上了哪个数字将这一列的非主元项目削减为0（在LU分解中用于L矩阵的产生）
        //FirstTransFormTimes表示做了几次第一类变换（行变换）

        FirstTransFormTimes = 0; //初始化为0
        size_t Row = mat.GetNumRow();
        size_t Col = mat.GetNumCol();

        size_t row_counter = 0;
        size_t col_counter = 0;
        T PivotElement;

        for (col_counter = 0; col_counter < Col; col_counter++)
        {
            row_counter = col_counter;


            //对角线主元是否为0？如果是则交换，如果不是则继续
            //如果需要交换，则采用初等变换进行

            PivotElement = mat(row_counter, row_counter);
            size_t tmp = 0;
            while (abs(PivotElement) < EPS)
            {
                // 如果主元为0则在该列一直寻找，直到找到一个主元不为0的数字
                tmp++;

                if (row_counter + tmp == Row)
                {
                    //表示本列中所有元素均为0
                    break;
                }
                PivotElement = mat(row_counter + tmp, row_counter);
            }

            if (tmp != 0 && row_counter + tmp < Row)
            {
                //该行的对角主元确实为0
                //且在该列中确实发现不为0的元素

                //进行交换
                mat.FirstTypeTransForm(row_counter + tmp, row_counter);
                FirstTransFormTimes++;
            }

            if (row_counter + tmp == Row)
            {
                continue;
            }

            PivotElement = mat(row_counter, row_counter);



            for (size_t tmp_row_counter = row_counter + 1; tmp_row_counter < Col; tmp_row_counter++)
            {
                //从对角的下一行开始
                T num = mat(tmp_row_counter, col_counter);
                double times = num / PivotElement;
                // 倍数

                mat.ThirdTypeTransForm(row_counter, tmp_row_counter, -times);
                RowReturn.push_back(tmp_row_counter);
                ColReturn.push_back(col_counter);
                NumReturn.push_back(-times);
            }
        }
        return mat;
    }


    static void givens(double a,double b,double *c,double *s)
    {
        //算法来源为 《数值线性代数》 第二版 P90
        if(abs(b) < EPS)
        {
            *c = 1;
            *s = 0;
            return;
        }
        else
        {
            double t;
            if(abs(b) > abs(a))
            {
                t = a / b;
                *s = 1 /sqrt(1+t*t);
                *c = (*s) *t;
                return;
            }
            else
            {
                t = b /a;
                *c = 1 /sqrt(1+t*t);
                *s = (*c) * t;
                return;
            }
        }
    }

    static const Matrix<T> LowtriSolve(const Matrix<T>& A_Input,const Matrix<T>& b_Input)
    {

        //解下三角形矩阵
        //采用前代法

        //理论参见 《数值线性代数》 第二版 P12
        if(!A_Input.isLowTri())
        {
            throw runtime_error("A is not LowTri");
        }

        Matrix<T> A(A_Input);
        Matrix<T> b(b_Input);


        Matrix<T> ans;

        size_t row = A.GetNumRow();

        ans.Resize(b.GetNumRow(), 1);
        for (size_t i = 0; i < A.GetNumRow() - 1; i++)
        {
            Matrix<T> bn(b.GetNumRow() - i - 1, 1);
            Matrix<T> An(b.GetNumRow() - i - 1, 1);
            bn = b.ExtractBlock(i + 1, 0, bn.GetNumRow(), bn.GetNumCol());
            An = A.ExtractBlock(i + 1, i, An.GetNumRow(), An.GetNumCol());

            ans(i,0) = b(i,0) / A(i, i);

            bn -= An*ans(i, 0);

            b.SetBlock(i + 1, 0, bn.GetNumRow(), bn.GetNumCol(), bn);
        }

        ans(row - 1, 0) = b(row - 1, 0) / A(row - 1, row - 1);

        return ans;
    }

    static const Matrix<T> UptriSolve(const Matrix<T>& A_Input,const Matrix<T>& b_Input)
    {
        //解上三角矩阵
        //采用回代法

        //理论参见 《数值线性代数》 第二版 P13

        if(!A_Input.isUpTri())
        {
            throw runtime_error("A is not UpTri");
        }

        Matrix<T> ans;

        Matrix<T> A(A_Input);
        Matrix<T> b(b_Input);

        ans.Resize(b.GetNumRow(),1);


        for(size_t i = A.GetNumRow() - 1; i >= 1; i--)
        {
            Matrix<T> bn(i,1);
            Matrix<T> An(i,1);

            bn = b.ExtractBlock(0,0,bn.GetNumRow(),bn.GetNumCol());
            An = A.ExtractBlock(0,i,An.GetNumRow(),An.GetNumCol());


            ans(i,0) = b(i,0) / A(i,i);
            bn -= An * ans(i,0);

            b.SetBlock(0,0,bn.GetNumRow(),bn.GetNumCol(),bn);
        }

        ans(0,0) = b(0,0) / A(0,0);

        return ans;
    }

    static int sgn(const T& value)
    {
        if(value > 0)
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }


};







#endif //CPP_MATRIX_UTILITY_H
