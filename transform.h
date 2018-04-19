#ifndef __TRANSFORM_H__
#define __TRANSFORM_H__
#include <vector>
#include "Matrix.h"


//2017.09.07 尚未测试
template<typename DataType>
Matrix<DataType> ToRowEchelonForm(const Matrix<DataType>& mat, std::vector<int>& RowReturn,
	std::vector<int>& ColReturn, std::vector<double>& NumReturn, int& FirstTransFormTimes)
{
	//采用消元将一个矩阵变为行最简矩阵
	//RowReturn，ColReturn，NumReturn返回的是消元时对哪一行哪一列加上了哪个数字将这一列的非主元项目削减为0（在LU分解中用于L矩阵的产生）
	//FirstTransFormTimes表示做了几次第一类变换（行变换）

	FirstTransFormTimes = 0; //初始化为0
	int Row = mat.GetNumRow();
	int Col = mat.GetNumCol();
	Matrix<DataType> matOutput;
	//输出行最简矩阵
	matOutput = mat;

	int row_counter = 0;
	int col_counter = 0;
	DataType PivotElement;

	for (col_counter = 0; col_counter < Col; col_counter++)
	{
		row_counter = col_counter;


		//对角线主元是否为0？如果是则交换，如果不是则继续
		//如果需要交换，则采用初等变换进行

		PivotElement = matOutput(row_counter, row_counter);
		int tmp = 0;
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



		for (int tmp_row_counter = row_counter + 1; tmp_row_counter < Col; tmp_row_counter++)
		{
			//从对角的下一行开始
			DataType num = matOutput(tmp_row_counter, col_counter);
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


template<typename DataType>
Matrix<DataType>& ToRowEchelonForm(Matrix<DataType>& mat, std::vector<int>& RowReturn,
	std::vector<int>& ColReturn, std::vector<double>& NumReturn, int& FirstTransFormTimes)
{
	//采用消元将一个矩阵变为行最简矩阵
	//RowReturn，ColReturn，NumReturn返回的是消元时对哪一行哪一列加上了哪个数字将这一列的非主元项目削减为0（在LU分解中用于L矩阵的产生）
	//FirstTransFormTimes表示做了几次第一类变换（行变换）

	FirstTransFormTimes = 0; //初始化为0
	int Row = mat.GetNumRow();
	int Col = mat.GetNumCol();

	int row_counter = 0;
	int col_counter = 0;
	DataType PivotElement;

	for (col_counter = 0; col_counter < Col; col_counter++)
	{
		row_counter = col_counter;


		//对角线主元是否为0？如果是则交换，如果不是则继续
		//如果需要交换，则采用初等变换进行

		PivotElement = mat(row_counter, row_counter);
		int tmp = 0;
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



		for (int tmp_row_counter = row_counter + 1; tmp_row_counter < Col; tmp_row_counter++)
		{
			//从对角的下一行开始
			DataType num = mat(tmp_row_counter, col_counter);
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





#endif // DECOMPOSITION_H
