#include <iostream>
#include <map>
 
int main()
{
    std::multimap<int, char> dict
    {
        {1, 'A'},
        {2, 'B'},
        {2, 'C'},
        {2, 'D'},
        {4, 'E'},
        {3, 'F'},
        {2, 'T'}
    };
 
    auto range = dict.equal_range(2);
 
    for (auto i = range.first; i != range.second; ++i)
        std::cout << i->first << ": " << i->second << '\n';
}