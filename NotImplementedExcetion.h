//
// Created by Yangwenlu on 2019/4/24.
//

#ifndef CPP_MATRIX_NOTIMPLEMENTEDEXCETION_H
#define CPP_MATRIX_NOTIMPLEMENTEDEXCETION_H

#include <string>

using namespace std;

class NotImplementedExcetion
{
public:
    NotImplementedExcetion();
    NotImplementedExcetion(const string& str);
    const string& whats() const throw();

private:
    string info;


};


#endif //CPP_MATRIX_NOTIMPLEMENTEDEXCETION_H
