// A "main" function to perform lot of insertions/updates/deletions.

#include <iostream>
#include <cstdlib>

#include "PriorityVector.hpp"

int main() {

    // Initializations
    PriorityVector tab(512);
    srand(314128);

    // Perform lot of tests
    for (size_t nrun = 0; nrun < 3; nrun++) {

        size_t count = 0;
        size_t insert = 0, remove = 0, update = 0;
        std::cout << "Starting run " << nrun << std::endl;

        for (size_t n = 0; n < 100000; n++) {
            int index = rand() % tab.size();

            // Set or remove a priority
            if (!tab.hasPriority(index)) {
                tab.setPriority(index, rand());
                count++;
                update++;
            } else if (rand() % 2) {
                tab.setPriority(index, rand());
                insert++;
            } else {
                tab.removePriority(index);
                count--;
                remove++;
            }

            // Check if the heap is correct
            if (count != tab.count() || tab.hasErrors()) {
                std::cout << "Invalid heap at iteration " << n << std::endl;
                throw "Error";
            }
        }

        // Remove all the indices
        if (count) {
            size_t value = tab[tab.getFirstIndex()];
            for (size_t n = count; n > 0; n--) {
                size_t newValue = tab[tab.removeFirstIndex()];
                if (value < newValue) {
                    std::cout << "Invalid order" << std::endl;
                    throw "Error";
                }
                value = newValue;
            }
        }

        // Check the size
        if (tab.count()) {
            std::cout << "Invalid order" << std::endl;
            throw "Error";
        }

        // Display results
        std::cout << "-> success, with "
            << insert << " insert, "
            << remove << " remove, "
            << update << " update, "
            << count  << " sorted" << std::endl;
    }

    return 0;
}
