#include "PriorityVector.hpp"


// Constructors

PriorityVector::PriorityVector(size_t length) :
    std::vector<double>(length, 0.0), heap(), inHeap(length, 0)
{
    heap.reserve(size() + rootIndex());
    heap.resize(rootIndex());
}

PriorityVector::PriorityVector(const std::vector<double>& priorities) :
    std::vector<double>(priorities), heap(), inHeap()
{
    heap.reserve(size() + rootIndex()); // better to declare its length
    heap.resize(rootIndex());
    inHeap.resize(size());
    for (size_t id = 0; id < size(); id++) {
        inHeap[id] = heap.size();
        heap.push_back(id);
        bubbleUp(inHeap[id]);
    }
}


// Add and remove functions

size_t PriorityVector::removeFirstIndex() {
    size_t id = getFirstIndex();
    removePriority(id);
    return id;
}

void PriorityVector::setPriority(size_t id, double value) {
    if (id >= size()) {
        throw std::out_of_range("Priority index is out of range.");
    }
    std::vector<double>::operator[](id) = value;
    if (hasPriority(id)) {
        bubbleUp(inHeap[id]);
        bubbleDown(inHeap[id]);
    } else {
        inHeap[id] = heap.size();
        heap.push_back(id);
        bubbleUp(inHeap[id]);
    }
}

void PriorityVector::removePriority(size_t id) {
    if (!hasPriority(id)) {
        throw std::out_of_range("This index has no priority.");
    }
    size_t k = inHeap[id];
    size_t kLast = heap.size() - 1;
    if (k == kLast) {
        heap.pop_back();
        inHeap[id] = 0;
    } else {
        swap(k, kLast);
        heap.pop_back();
        inHeap[id] = 0;
        bubbleUp(k);
        bubbleDown(k);
    }
}

void PriorityVector::removePriority() {
    heap.resize(rootIndex());
    inHeap.assign(size(), 0);
}


// Debug function

size_t PriorityVector::hasErrors() const {
    size_t errors = 0;
    for (size_t k = rootIndex() + 1; k < heap.size(); k++) {
        if (isLessThan(parentIndex(k), k)) {
            errors++;
        }
    }
    for (size_t id = 0; id < size(); id++) {
        if (hasPriority(id) && heap[inHeap[id]] != id) {
            errors++;
        }
    }
    return errors;
}


// Swap and bubble functions

void PriorityVector::swap(size_t k, size_t other) {
    size_t id = heap[k];
    inHeap[heap[k] = heap[other]] = k;
    inHeap[heap[other] = id] = other;
}

void PriorityVector::bubbleUp(size_t k) {
    while (hasParent(k) && isLessThan(parentIndex(k), k)) {
        swap(k, parentIndex(k));
        k = parentIndex(k);
    }
}

void PriorityVector::bubbleDown(size_t k) {
    size_t child;
    while ((child = firstChildIndex(k)) < heap.size()) {
        if (child + 1 < heap.size() && isLessThan(child, child + 1)) {
            child = child + 1;
        }
        if (isLessThan(k, child)) {
            swap(k, child);
        } else {
            break;
        }
        k = child;
    }
}
