#include <iostream>
#include <vector>
using namespace std;
int partitation(vector<int> arr, int l, int r){
    int p = arr[r];
    int start = l, end = r-1;
    while(start != arr.size()){
        if(arr[start] > p) continue;
        start ++;
    }
    while(end != l){
        if(arr[end] < p) continue;
        end --;
    }
    swap(arr[start],arr[end]);
    if(end >= start)
            swap(arr[start],arr[end]);
}
void quickSort(vector<int> arr,int l, int r){
    if(l < r){
        int p = partitation(arr,l,r);
        quickSort(arr,l,p - 1);
        quickSort(arr,p,r);
    }
}

int main() {
	//code
	vector<int> arr;
    arr.push_back(23);
    arr.push_back(45);
    arr.push_back(2);
    arr.push_back(5);
    quickSort(arr,0,arr.size());
    cout<<arr[2]<<endl;
    return 0;
}
