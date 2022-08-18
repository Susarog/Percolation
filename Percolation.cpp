#include <iostream>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <array>
#include <stdlib.h>
#include <time.h>
using namespace std;

class UnionFind {
public:
    UnionFind(int n = 0) {
        for (int i = 0; i < n; i++) {
            tree.push_back(-1);
        }
    }
    void validate(int v1) {
        if (v1 > tree.size() || v1 < 0) {
            throw std::invalid_argument("invalid index");
        }
    }

    int sizeOf(int v1) {
        validate(v1);
        return -(tree.at(find(v1)));
    }
    int parent(int v1) {
        validate(v1);
        return tree.at(v1);
    }
    bool connected(int v1, int v2) {
        return find(v1) == find(v2);
    }

    void union1(int v1, int v2) {
        if (connected(v1, v2)) {
            return;
        }
        int weight1 = sizeOf(v1);
        int weight2 = sizeOf(v2);
        int compare = weight1 - weight2;
        //v1 has more weight
        if (compare > 0) {
            tree.at(find(v2)) = find(v1);
            tree.at(find(v1)) = -(weight1 + weight2);
        }
        //v2 has more weight
        else if (compare < 0) {
            tree.at(find(v1)) = find(v2);
            tree.at(find(v2)) = -(weight1 + weight2);
        }
        //same weight
        else {
            tree.at(find(v2)) = find(v1);
            tree.at(find(v1)) = -(weight1 + weight2);
        }
    }
    //returns the index of the root
    int find(int v1) {
        validate(v1);
        while (true) {
            if (tree.at(v1) < 0) {
                break;
            }
            v1 = tree.at(v1);
        }
        return v1;
    }
    int size() {
        return tree.size();
    }
    void print() {
        for (int i = 0; i < tree.size(); i++) {
            std::cout << tree.at(i) << std::endl;
        }
    }
private:
    std::vector<int> tree;
};

class Percolation {
public:
    Percolation(int N = 10) { // create N-by-N grid, with all sites initially blocked
        if (N <= 0) {
            throw std::invalid_argument("Input is less than or equal to 0");
        }
        disjoint = UnionFind(N * N + 2);
        size = N;
        //create a 2d array with 2 secret slots that the user cannot get to.

        for (int i = 0; i < N; i++) {
            vector<int> temp;
            for (int j = 0; j < N; j++) {
                temp.push_back(0);
            }
            percolation.push_back(temp);
        }
        vector<int> temp;
        topVirtualVal = N * N;
        bottomVirtualVal = N * N + 1;
        temp.push_back(topVirtualVal);
        temp.push_back(bottomVirtualVal);
        percolation.push_back(temp);
        //union upper row
        for (int i = 0; i < N; i++) {
            disjoint.union1(i, topVirtualVal);
        }
        //union bottom row
        for (int i = N * N - N; i < N * N; i++) {
            disjoint.union1(i, bottomVirtualVal);
        }

    }
    // open the site (row, col) if it is not open already
    //check for neighbors on all 4 directions
    void open(int row, int col) {
        if (this->isOpen(row, col)) {
            return;
        }
        counter++;
        this->percolation[row][col] = 1;
        int uniqueVal = uniqueValue(row, col);

        // check for all edges of opened percolation check if row or col becomes >= size then check if there is an open slot
        if (row + 1 < size && isOpen(row + 1, col)) {
            int secondVal = this->uniqueValue(row + 1, col);
            disjoint.union1(uniqueVal, secondVal);
        }
        if (row - 1 >= 0 && isOpen(row - 1, col)) {
            int secondVal = this->uniqueValue(row - 1, col);
            disjoint.union1(uniqueVal, secondVal);
        }
        if (col + 1 < size && isOpen(row, col + 1)) {
            int secondVal = this->uniqueValue(row, col + 1);
            disjoint.union1(uniqueVal, secondVal);
        }
        if (col - 1 >= 0 && isOpen(row, col - 1)) {
            int secondVal = this->uniqueValue(row, col - 1);
            disjoint.union1(uniqueVal, secondVal);
        }
        if (isFull(row, col)) {
            percolationBruh(row, col);
        }

    }
    bool isOpen(int row, int col) {  // is the site (row, col) open?
        isInvalid(row, col);
        return (percolation[row][col] > 0);
    }
    bool isFull(int row, int col) {  // is the site (row, col) full?
        isInvalid(row, col);
        int val = this->uniqueValue(row, col);
        return disjoint.connected(val, topVirtualVal);
    }
    int numberOfOpenSites() {        // number of open sites
        return counter;
    }
    bool percolates() {              // does the system percolate? (check top and bottom of 2d array)
        return disjoint.connected(topVirtualVal, bottomVirtualVal);
    }
    //worst case can be around n^2 need to find a better way to change the value to 2
    void percolationBruh(int row, int col) {
        this->percolation[row][col] = 2;
        if (row + 1 < size && isOpen(row + 1, col) && this->percolation[row + 1][col] != 2) {
            percolationBruh(row + 1, col);
        }
        if (row - 1 >= 0 && isOpen(row - 1, col) && this->percolation[row - 1][col] != 2) {
            percolationBruh(row - 1, col);
        }
        if (col + 1 < size && isOpen(row, col + 1) && this->percolation[row][col + 1] != 2) {
            percolationBruh(row, col + 1);
        }
        if (col - 1 >= 0 && isOpen(row, col - 1) && this->percolation[row][col - 1] != 2) {
            percolationBruh(row, col - 1);

        }
    }
    void print() {
        for (int i = 0; i < size; i++) {
            cout << "{ ";
            for (int j = 0; j < size; j++) {
                cout << percolation[i][j];
            }
            cout << " }\n";
        }
        cout << endl;
    }
    void printUnionFind() {
        disjoint.print();
    }
private:
    vector<vector<int>> percolation;
    UnionFind disjoint;
    int size = 0;
    int counter = 0;
    int topVirtualVal = 0;
    int bottomVirtualVal = 0;
    int uniqueValue(int i, int j) {
        return size * i + j;
    }
    void isInvalid(int row, int col) {
        if (row > size || row < 0 || col > size || col < 0) {
            throw std::invalid_argument("IndexOutOfBoundsException");
        }
    }
};
class PercolationFactory {
public:
    Percolation make(int N) {
        return Percolation(N);
    }
};
class PercolationStats {
public:
    // perform T independent experiments on an N-by-N grid
    PercolationStats(int N, int T, PercolationFactory Pf) {
        size = N;
        numsOfExp = T;
        bruh = Pf;
        meanVal = 0;
        stdDevVal = 0;
    }
    void simulation1() {
        vector<double> tracker(0);
        double sum = 0;
        int randRow = 0;
        int randCol = 0;
        srand(time(NULL));
        Percolation idk;
        for (int i = 0; i < numsOfExp; i++) {
            idk = bruh.make(size);
            while (!idk.percolates()) {
                randRow = rand() % size;
                randCol = rand() % size;
                if (!idk.isOpen(randRow, randCol)) {
                    idk.open(randRow, randCol);
                }
            }
            tracker.push_back((double)idk.numberOfOpenSites() / ((double)size * size));
            sum += (double)idk.numberOfOpenSites() / ((double)size * size);
        }
        meanVal = (double)sum / (double)numsOfExp;
        sum = 0;
        for (int i = 0; i < numsOfExp; i++) {
            sum += pow((tracker.at(i) - meanVal), 2);
        }
        stdDevVal = sqrt((sum / (numsOfExp - 1.0)));
    }
    // sample mean of percolation threshold
    double mean() {
        return meanVal;
    }
    // sample standard deviation of percolation threshold
    double stddev() {
        return stdDevVal;
    }
    // low endpoint of 95% confidence interval

    double confidenceLow() {
        int mean = this->mean();
        return mean - (1.96 * this->stddev()) / sqrt(numsOfExp);
    }
    // high endpoint of 95% confidence interval

    double confidenceHigh() {
        int mean = this->mean();
        return mean + (1.96 * this->stddev()) / sqrt(numsOfExp);
    }
private:
    PercolationFactory bruh;
    int numsOfExp;
    int size;
    double meanVal;
    double stdDevVal;
};


int main()
{
    PercolationFactory percFactory;
    PercolationStats percStats(20, 400, percFactory);
    percStats.simulation1();
    cout << "mean: " << percStats.mean() << endl << "stdDeviation: " << percStats.stddev() << endl << "confidenceHigh: " << percStats.confidenceHigh() << endl << "confidence Low: " << percStats.confidenceLow();
}
