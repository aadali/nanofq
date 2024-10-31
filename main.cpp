#include "test.h"
#include "Adapter.h"

using namespace std;

int main()
{
    cout << "hello world" << endl;
    auto a = barcode_info::get_trim_info();
    SequenceInfo& seq_info  {a.find("SQK-LSK114")->second};
    cout << "hello world" << endl;
}
