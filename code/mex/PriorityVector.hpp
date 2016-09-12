/* PriorityVector
 * This class implements a vector of priorities.
 *
 * It allows:
 * - the storage of a vector of priorities.
 * - an instant access to highest priority index.
 * - a fast priorities update.
 * - a fast priorities removal.
 *
 * Guillaume Tartavel, June 2012
 */


#ifndef _PRIORITY_VECTOR_HPP_
#define _PRIORITY_VECTOR_HPP_


#include <vector>
#include <stdexcept>


class PriorityVector : private std::vector<double> {

    // Arguments
    std::vector<size_t> heap;
    std::vector<size_t> inHeap;

public:

    // Constructors
    PriorityVector(size_t length);
    PriorityVector(const std::vector<double>& priorities);

    // Sizes
    inline size_t size() const
        { return std::vector<double>::size(); }
    inline size_t count() const
        { return heap.size() - rootIndex(); }

    // Index of the highest priorities
    inline size_t getFirstIndex() const
        { return heap[rootIndex()]; }
    size_t removeFirstIndex();

    // Get and set priorities
    inline double operator[](size_t id) const
        { return std::vector<double>::operator[](id); }
    inline double getPriority(size_t id) const
        { return operator[](id); }
    void setPriority(size_t id, double value);

    // Manage undefined priorities
    void removePriority(size_t id);
    void removePriority();
    inline bool hasPriority(size_t id) const
        { return (inHeap[id] > 0); }
    inline const std::vector<size_t>& hasPriority() const
        { return inHeap; }

    // Debug function
    size_t hasErrors() const;

private:

    // Swap and bubble functions
    void swap(size_t k, size_t other);
    void bubbleUp(size_t id);
    void bubbleDown(size_t id);

    // Heap relationships
    inline bool hasParent(size_t k) const
        { return (k > rootIndex()); }
    inline size_t rootIndex() const
        { return 1; }
    inline size_t parentIndex(size_t k) const
        { return k / 2; }
    inline size_t firstChildIndex(size_t k) const
        { return k * 2; }
    inline bool isLessThan(size_t k, size_t other) const
        { return (getPriority(heap[k]) < getPriority(heap[other])); }
};


#endif /* _PRIORITY_VECTOR_HPP_ */
