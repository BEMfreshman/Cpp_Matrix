//
// Created by Yangwenlu on 2019/4/24.
//

#include "NotImplementedExcetion.h"


NotImplementedExcetion::NotImplementedExcetion()
{

}

NotImplementedExcetion::NotImplementedExcetion(const string &str)
:info(str)
{

}

const string& NotImplementedExcetion::whats() const throw()
{
    return info;
}