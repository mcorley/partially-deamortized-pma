#include <iostream>
#include <stdint.h>
#include <cmath>
#include <vector>
#include "pma.h"
using namespace std;

static void dump_stats(const pma& database)
{
  cout << "Capacity: "  << database.capacity()            << endl;
  cout << "Size: "      << database.size()                << endl;
  cout << "SegSize: "   << database.segment_size()        << endl;
  cout << "Segments: "  << database.number_of_segments()  << endl;
  cout << "Height: "    << database.tree_height()         << endl;
}

static void dump_pma_contents(const pma& database)
{
  cout << "pma:  [";
  for (uint32_t i = 0; i < database.capacity(); ++i)
    cout << database[i] << "|";
  cout << "]" << endl;
}

static void dump_pma_free_index_bitmap(const pma& database)
{
  cout << "free: [";
  for (uint32_t i = 0; i < database.capacity(); ++i) {
    if (database.index_is_free(i)) cout << "0" << "|";
    else cout << "1" << "|";
  }
  cout << "]" << endl;
}

static void dump_upper_density_thresholds(const pma& database)
{
  cout << "UDTs : \n";
  for (int i = 0; i <= database.tree_height(); ++i)
    cout << database.upper_density_threshold(i) << endl;
  cout << endl;
}

int main(int argc, char *argv[]) 
{
  pma database;
  dump_stats(database);
  dump_pma_contents(database);
  dump_pma_free_index_bitmap(database);
  cout << endl;
  for (uint32_t i = 0; i < 4; ++i) {
    database.insert(i);
    dump_pma_contents(database);
    dump_pma_free_index_bitmap(database);
    cout << endl;
  }
  return 0;
}
