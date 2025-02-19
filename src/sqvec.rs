use std::alloc::{self, Layout, alloc, dealloc, realloc};
use std::fmt::Debug;
use std::marker::PhantomData;
use std::ops::{Index, IndexMut};
use std::ptr::NonNull;

// todo
// - add shrinking the array to size
// - efficient item drop
// - c interop
// - size checks/debug assertions
// - unique or phantomdata?
// - initialisation macro
// - clone trait
// - from slice/box trait
// - intoiter
// - insert
// - remove
// - store seglen instead of log_seglen?
// - zero sized types
// - dropck/phantomdata?
// - benchmark different iteration strategies
// - benchmark random operations against vec
// - iter_mut()
// - usize as index?

/// A non-contiguos extensible array that achieves asymptotically optimal
/// space overhead of O(√n)). Specification: <https://jmmackenzie.io/pdf/mm22-adcs.pdf>
pub struct SqVec<T> {
    dope: NonNull<NonNull<T>>,
    dope_cap: u32,
    len: u32,
    alloc_seg_count: u32,
    log_seglen: u32,
    _owns_t: PhantomData<T>,
}

impl<T: Debug> Debug for SqVec<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl<T> SqVec<T> {
    /// Constructs a new, empty `SqVec<T>`.
    ///
    /// The sqvec will not allocate until elements are pushed onto it.
    pub fn new() -> Self {
        assert!(
            std::mem::size_of::<T>() != 0,
            "This data structure does not support ZST's"
        );
        Self {
            dope: NonNull::dangling(),
            dope_cap: 0,
            len: 0,
            alloc_seg_count: 0,
            log_seglen: 1,
            _owns_t: PhantomData,
        }
    }

    /// grows the dope vector to the correct size to fit the next set of segments
    fn grow_dope(&mut self) {
        let (new_cap, new_layout) = if self.dope_cap == 0 {
            (2, unsafe {
                Layout::array::<NonNull<T>>(2).unwrap_unchecked()
            })
        } else {
            let new_cap = self.dope_cap + Self::segs_of_len(self.log_seglen);
            let new_layout =
                Layout::array::<NonNull<T>>(new_cap as usize).expect("Allocation too large");
            (new_cap, new_layout)
        };

        let handle = if self.dope_cap == 0 {
            unsafe { alloc(new_layout) }
        } else {
            let old_layout =
                unsafe { Layout::array::<NonNull<T>>(self.dope_cap as usize).unwrap_unchecked() };
            let old_handle = self.dope.as_ptr() as *mut u8;
            unsafe { realloc(old_handle, old_layout, new_layout.size()) }
        };

        self.dope = NonNull::new(handle as *mut NonNull<T>)
            .unwrap_or_else(|| alloc::handle_alloc_error(new_layout));
        self.dope_cap = new_cap;
    }

    /// Allocates a new segment of length `1 << self.log_seglen` at index `segnum` in dope.
    ///
    /// This function assumes no segment is allocated at `segnum`.
    fn alloc_seg(&mut self) {
        let layout = Layout::array::<T>(1 << self.log_seglen).expect("Allocation too large");
        let handle = unsafe { alloc(layout) };
        let ptr =
            NonNull::new(handle as *mut T).unwrap_or_else(|| alloc::handle_alloc_error(layout));
        unsafe { self.dope.add(self.alloc_seg_count as usize).write(ptr) };
        self.alloc_seg_count += 1;
    }

    /// Adds `val` to the end of the SqVec, increasing its length by 1.
    pub fn push(&mut self, val: T) {
        let (segnum, offset) = Self::mapping(self.len);
        if offset == 0 {
            if segnum == Self::last_seg(self.log_seglen) + 1 {
                self.log_seglen += 1;
            }
            if segnum == self.alloc_seg_count {
                if segnum == self.dope_cap {
                    self.grow_dope();
                }
                self.alloc_seg();
            }
        }
        unsafe { self.item_ptr_inner(segnum, offset).write(val) };
        self.len += 1;
    }

    /// Removes the last value on the SqVec and returns it in an Option.
    /// If the SqVec is empty, Returns None
    pub fn pop(&mut self) -> Option<T> {
        if self.len == 0 {
            return None;
        }
        self.len -= 1;
        let (segnum, offset) = Self::mapping(self.len);
        if (segnum == Self::first_seg(self.log_seglen)) && (offset == 0) {
            self.log_seglen = 1.max(self.log_seglen - 1);
        }
        Some(unsafe { self.item_ptr_inner(segnum, offset).read() })
    }

    /// swaps the elements at index's ix1 and ix2
    pub fn swap(&mut self, ix1: u32, ix2: u32) {
        assert!(ix1 < self.len && ix2 < self.len, "Index out of range");
        unsafe { self.swap_unchecked(ix1, ix2) };
    }

    pub unsafe fn swap_unchecked(&mut self, ix1: u32, ix2: u32) {
        unsafe { std::mem::swap(self.item_ptr(ix1).as_mut(), self.item_ptr(ix2).as_mut()) };
    }

    /// Returns the length of the SqVec (number of elements)
    #[inline(always)]
    pub fn len(&self) -> u32 {
        self.len
    }

    /// Returns true if the sqvec is empty and false otherwise
    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Maps the given index to the segment number and index within the
    /// segment (offset) where the element at index `v` can be found.
    #[inline(always)]
    fn mapping(v: u32) -> (u32, u32) {
        let b = if v == 0 {
            1
        } else {
            (33 - v.leading_zeros()) >> 1
        };
        let segnum = (v >> b) + (1 << (b - 1)) - 1;
        let offset = v & ((1 << b) - 1);
        (segnum, offset)
    }

    /// Gives the dope index of the first segment of length (1 << log_seglen)
    #[inline(always)]
    fn first_seg(log_seglen: u32) -> u32 {
        if std::intrinsics::unlikely(log_seglen == 1) {
            0
        } else {
            3 * (1 << (log_seglen - 2)) - 1
        }
    }

    /// Gives the dope index of the last segment of length (1 << log_seglen)
    #[inline(always)]
    fn last_seg(log_seglen: u32) -> u32 {
        3 * (1 << (log_seglen - 1)) - 2
    }

    /// Returns the max number of segments of lenght (1 << log_seglen) that can
    /// be allocated
    #[inline(always)]
    fn segs_of_len(log_seglen: u32) -> u32 {
        if std::intrinsics::unlikely(log_seglen == 1) {
            2
        } else {
            3 * (1 << (log_seglen - 2))
        }
    }

    /// Returns the pointer to the item in the SqVec at index `ix` (Does not check if ix is < self.len)
    #[inline(always)]
    unsafe fn item_ptr(&self, ix: u32) -> NonNull<T> {
        let (segnum, offset) = Self::mapping(ix);
        unsafe { self.item_ptr_inner(segnum, offset) }
    }

    ///
    #[inline(always)]
    unsafe fn item_ptr_inner(&self, segnum: u32, offset: u32) -> NonNull<T> {
        unsafe { self.dope.add(segnum as usize).read().add(offset as usize) }
    }

    fn iter(&self) -> Iter<T> {
        Iter {
            data: &self,
            seglen: 2,
            segnum: 0,
            offset: 0,
            last_seg_ix_of_len: 1,
            ix: 0,
        }
    }

    fn t_iter(&self) -> TIter<T> {
        TIter { data: &self, ix: 0 }
    }
}

pub struct IterMut<'a, T> {
    data: &'a SqVec<T>,
    ix: u32,
}

impl<'a, T> Iterator for IterMut<'a, T> {
    type Item = &'a mut T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.ix == self.data.len() {
            return None;
        }
        let out = unsafe { self.data.item_ptr(self.ix).as_mut() };
        self.ix += 1;
        return Some(out);
    }
}

pub struct TIter<'a, T> {
    data: &'a SqVec<T>,
    ix: u32,
}

impl<'a, T> Iterator for TIter<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.ix == self.data.len() {
            return None;
        }
        let out = unsafe { self.data.item_ptr(self.ix).as_ref() };
        self.ix += 1;
        return Some(out);
    }
}

pub struct Iter<'a, T> {
    data: &'a SqVec<T>,
    seglen: u32,
    segnum: u32,
    offset: u32,
    last_seg_ix_of_len: u32,
    ix: u32,
}

impl<'a, T> Iterator for Iter<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.ix == self.data.len() {
            return None;
        }

        if self.offset == self.seglen {
            // next segnum... also check wether seglen should double
            if self.segnum == self.last_seg_ix_of_len {
                self.last_seg_ix_of_len = 3 * self.seglen - 2;
                self.seglen <<= 1;
            }
            self.segnum += 1;
            self.offset = 0;
        }

        let out = unsafe { self.data.item_ptr_inner(self.segnum, self.offset).as_ref() };
        self.offset += 1;
        self.ix += 1;
        return Some(out);
    }
}

impl<T> Index<u32> for SqVec<T> {
    type Output = T;

    fn index(&self, ix: u32) -> &Self::Output {
        assert!(ix < self.len, "Index out of bounds.");
        unsafe { self.item_ptr(ix).as_ref() }
    }
}

impl<T> IndexMut<u32> for SqVec<T> {
    fn index_mut(&mut self, ix: u32) -> &mut Self::Output {
        assert!(ix < self.len, "Index out of bounds.");
        unsafe { self.item_ptr(ix).as_mut() }
    }
}

impl<T> Drop for SqVec<T> {
    fn drop(&mut self) {
        if self.dope_cap == 0 {
            return;
        }
        // drop items in SqVec
        for i in 0..self.len {
            unsafe { self.item_ptr(i).drop_in_place() };
        }

        let dope_slice = unsafe {
            std::slice::from_raw_parts(self.dope.as_ptr(), self.alloc_seg_count as usize)
        };
        let mut ix1 = 0_usize;
        let max_lsl = u32::ilog2((self.dope_cap + 2) / 3) + 1;
        for lsl in 1..=max_lsl {
            let ix2 = Self::last_seg(lsl).min(self.alloc_seg_count.saturating_sub(1)) as usize;
            let layout = Layout::array::<T>(1 << lsl).unwrap();
            for &ptr in dope_slice[ix1..=ix2].iter() {
                unsafe { dealloc(ptr.as_ptr() as *mut u8, layout) };
            }
            ix1 = ix2 + 1;
        }

        unsafe {
            dealloc(
                self.dope.as_ptr() as *mut u8,
                Layout::array::<NonNull<T>>(self.dope_cap as usize).unwrap(),
            )
        };
    }
}

unsafe impl<T: Send> Send for SqVec<T> {}
unsafe impl<T: Sync> Sync for SqVec<T> {}

#[cfg(test)]
mod tests {
    use super::SqVec;
    use rand::{Rng, random, thread_rng};

    #[test]
    fn single_box_drop() {
        let mut sqvec = SqVec::<Box<u32>>::new();
        sqvec.push(Box::new(32));
        sqvec.pop();
    }

    #[test]
    fn drop_len_lt_cap() {
        let mut sqvec = SqVec::<Box<u32>>::new();
        for i in 0..200 {
            sqvec.push(Box::new(i));
        }
        for _ in 0..170 {
            sqvec.pop();
        }
    }

    #[test]
    fn drop_empty() {
        let _sqvec = SqVec::<Box<u32>>::new();
    }

    #[test]
    /// randomised test of mutating api
    fn fuzz() {
        let mut rng = thread_rng();
        let mut vec = Vec::<Box<u32>>::new();
        let mut sqvec = SqVec::<Box<u32>>::new();
        for _ in 0u32..2048 {
            let choice: u8 = rng.gen_range(0..=4);
            match choice {
                0 => {
                    if !vec.is_empty() {
                        let ix = rng.gen_range(0..vec.len());
                        assert_eq!(vec[ix], sqvec[ix as u32], "Not equal at index: {ix}");
                    }
                }
                1 => {
                    assert_eq!(vec.pop(), sqvec.pop());
                }
                2..=4 => {
                    let val: Box<u32> = Box::new(random());
                    vec.push(val.clone());
                    sqvec.push(val.clone());
                }
                _ => (),
            };
            assert_eq!(vec.len(), sqvec.len() as usize);
        }
    }

    #[test]
    fn lifetime() {
        let mut sqvec = SqVec::<Box<u32>>::new();
        sqvec.push(Box::new(23));
        let a = sqvec.pop().unwrap();
        drop(a);
        sqvec.push(Box::new(324));
    }

    #[test]
    fn iteration() {
        let mut rng = thread_rng();
        let mut sqvec = SqVec::<u32>::new();
        let mut vec = Vec::<u32>::new();
        for _ in 0u32..2048 {
            let a = rng.r#gen();
            sqvec.push(a);
            vec.push(a);
        }
        let mut sqiter = sqvec.iter();
        let mut viter = vec.iter();
        for i in 0u32..2048 {
            assert_eq!(sqiter.next(), viter.next(), "{i}");
        }
    }
}
