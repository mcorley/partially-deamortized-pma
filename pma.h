
#ifndef PMA_H
#define PMA_H

/** 
 * PMA  Packed-Memory Array
 * A cache-oblivious solution for maintaining a dynamic set of elements in 
 * sorted order. It supports operations insert, delete, and scan. Because 
 * the elements are stored physically in sorted order in memory or on disk, 
 * the PMA can be used to support extremly efficient range queries.
 *
 * The PMA is divided into $/Theta(N/logn)$ segments, and a contiguous run of
 * segments is called a window. The PMA is viewed in terms of a tree structure, 
 * where nodes of the tree are windows. The root node is the window containing 
 * all segments and a leaf node is a window containing a single segment. The
 * tree is implicitly rather than explicitly maintained.
 */
class pma {
  public:
    // The initial number of elements in the packed-memory array to maintain.
    // Should be a power of two.
    static const int INITIAL_CAPACITY = 4;
    static const int SCALE_FACTOR = 2;

    // Constant minimium and maximum densities.
    //
    // THRESHOLDS
    // The nodes at each height h have an upper density threshold t_h and a
    // lower density threshold p_h, which together determine the acceptable
    // density of keys within a window of 2^h segments. As node height
    // increases the udts decrease and ldts increase.
    //           D_min = p_0 <...< p_h < t_h <...< t_0 = D_max
    static const double LEAF_LOWER_DENSITY = 0.1;
    static const double ROOT_LOWER_DENSITY = 0.2;
    static const double ROOT_UPPER_DENSITY = 0.5;
    static const double LEAF_UPPER_DENSITY = 1.0;

  private:
    // The height of the root i.e the height of the tree.  
    int _implicit_tree_height;
    
    // The number of array positions in each leaf node i.e. segment.
    uint32_t _segment_size;

    // The number of elements in the packed-memory array.
    uint32_t _size;

    // A given bit in the bitmap is set if the corresponding index in the pma
    // is free, and clear if the corresponding index in the pma is in use.
    std::vector<bool> _free_index_bitmap;

    // The allocated storage space for the elements of the pma.
    std::vector<int> _storage;

  public:
    /** 
     * Default constructor: 
     * Constructs an empty packed-memory array, with no content and a 
     * size of zero. 
     */
    pma();

    /**
     * Destructs the packed-memory array. This calls each of the contained 
     * element's destructors, and deallocates all the storage capacity allocated 
     * by the vector.
     */
    ~pma();

    /**
     * Returns a reference to the element at position n in the packed-memory 
     * array.
     */
    int& operator[] (uint32_t& n);
    const int& operator[] (uint32_t& n) const;

    /**
     * Returns the size of the allocated storage space for the elements of the 
     * packed-memory array. The capacity is not necessarily equal to the number 
     * of elements that conform the underlying pma content (this can be obtained 
     * with member pma::size), but the capacity of the actual allocated space, 
     * which is either equal or greater than the content size.
     */
    uint32_t capacity() const;

    /**
     * Scans the window of the packed-memory array starting at index window
     * through index window+length and clears the storage contents and free
     * index bitmap information essentially nullifying this window.
     */
    void clear_window(const uint32_t& window, const uint32_t& length);

    /**
     * Removes from the packed-memory array a single element (x). This 
     * effectively reduces the pma size by the number of elements removed, 
     * calling each element's destructor before.
     */
    void erase(const int& x);

    /**
     * Returns whether index at position indexno in the free_index_bitmap is set.
     */
    bool index_is_free(uint32_t indexno) const;

    /**
     * The packed-memory array is extended by inserting new elements. This 
     * effectively increases the pma size, which causes an automatic reallocation 
     * of the allocated storage space if, and only if, the new pma size surpasses 
     * the current pma ROOT_UPPER_DENSITY. Rebalances of the pma may also be
     * triggered as a result of an insertion.
     */
    void insert(const int& x);

    /**
     * Computes the lower density threshold for a window at a given height in
     * the tree. As node height increases, the lower density threshold
     * increases. The threshold for nodes at height l are defined as:
     *                  p_l = p_h - (p_h - p_0)(h - l)/h
     * where h is the height of the tree, and 0 is at a leaf node.
     */
    double lower_density_threshold(int height) const;

    /**
     * Returns the number of segments i.e. leaf nodes in the tree. The number
     * of segments is always a power of two.
     */
    int number_of_segments() const;

    /**
     * Returns the index in the given segment of the packed-memory array to 
     * insert x into.
     * @param segment The index that starts the segment
     * @param x       The value of the element to be inserted.
     */
    uint32_t position_to_insert(const uint32_t& segment, const int& x) const;

    /**
     * Returns the index in the packed-memory array that holds the immediate
     * predecessor of x. Returns 0 if the pma is empty.
     */
    uint32_t predecessor(const int& x) const;

    /**
     * Returns the index in the packed-memory array that starts the segment
     * (leaf node) to insert x into.
     */
    uint32_t segment_to_insert(const int& x) const;    

    /**
     * We rebalance a node u_h of height h if u_h is within threshold, but we
     * detect that a child node u_h-1 is outside of threshold. Rebalances are
     * triggered by inserts or deletes that push one descendent node at each
     * height above its upper threshold t_i or below its lower threshold p_i.
     * @param segment The index that starts the segment out of balance.
     */
    void rebalance(const uint32_t& segment);
    void naive_rebalance(const uint32_t& window, const uint32_t& length);    
    void one_phase_rebalance(const uint32_t& window, const uint32_t& length);
    
    /**
     * When the packed-memory array becomes too full or too empty we recopy the
     * elements into a new pma that is a constant factor larger or smaller.
     */
    void resize();

    /**
     * Returns the number of indexes that are contained within a single 
     * segment.
     */
    uint32_t segment_size() const;

    /**
     * Returns the number of elements that conform the pma's content. This is 
     * the number of actual objects held in the pma, which is not necessarily 
     * equal to its storage capacity. 
     */
    uint32_t size() const;

    /**
     * Returns the height of the tree i.e how many implicit levels we are
     * maintaing. 
     */
    int tree_height() const;

    /**
     * Computes the lower density threshold for a window at a given height in
     * the tree. As node height increases, the lower density threshold
     * increases. The threshold for nodes at height l are defined as:
     *                  p_l = p_h - (p_h - p_0)(h - l)/h
     * where h is the height of the tree, and 0 is at leaf node level.
     */
    double upper_density_threshold(int height) const;

    /**
     * Returns the number of array positions in a node of height h. For leaf
     * nodes, i.e. segments of size one, the capacity is equal to the segment
     * size. For each level up in the tree, because we maintain segment sizes
     * that are powers of two, we can do a simple bit shift of the segment size
     * with the height of the tree we wish to compute for. Leaf nodes are of
     * height 0 and the root is of height _implicit_tree_height.
     *
     *  e.g. if we have segments of size 4 with an implicit tree height of 4,
     *    (leaf) cap at height 0 : 4 << 0 = 4
     *           cap at height 1 : 4 << 1 = 8
     *           cap at height 2 : 4 << 2 = 16
     *    (root) cap at height 3 : 4 << 3 = 32
     */
    uint32_t window_capacity(int height) const;

    /**
     * Scans the window of the packed-memory array starting at index window
     * through index window+length and returns the number of elements contained.
     */
    uint32_t window_size(const uint32_t& window, const uint32_t length) const;

    /**
     * Returns whether a node of height h has children that are inside their
     * density thresholds.
     */
    bool within_balance() const;    

  private:
    /**
     * Computes the next highest power of 2 of 32-bit value v.
     * From Bit Twiddling Hacks by Sean Eron Anderson
     * http://graphics.stanford.edu/~seander/bithacks.html
     */
    static inline uint32_t next_power_of_2(uint32_t v)
    {
      v--;
      v |= v >> 1;
      v |= v >> 2;
      v |= v >> 4;
      v |= v >> 8;
      v |= v >> 16;
      v++;
      return v;
    }

    /** Is v a power of 2? */
    static inline bool is_power_of_2(uint32_t v) {
      return v & (v - 1);
    }
};

#endif // PMA_H
