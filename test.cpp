#include <iostream>

using namespace std;
int main() {

    int f, fcmax;
    fcmax = 6;
    f = fcmax % 2 == 0 ? fcmax / 2 : fcmax / 2 + 1;
    cout << f << endl;

    fcmax = 7;
    f = fcmax % 2 == 0 ? fcmax / 2 : fcmax / 2 + 1;
    cout << f << endl;

    return 0;
}