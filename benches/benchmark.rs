use criterion::{Criterion, black_box, criterion_group, criterion_main};
use rand::{Rng, rng};
use sqvec::SqVec;

fn push_benchmark_sqvec(c: &mut Criterion) {
    // let mut rng = rng();

    c.bench_function("sqvec push", |c| {
        c.iter(|| {
            let mut sqvec = SqVec::<u32>::new();
            for i in 0..(1 << 23) {
                sqvec.push(i);
            }
        });
    });
}

fn sqvec_iter_bench1(c: &mut Criterion) {
    let mut sqvec = SqVec::<u32>::new();
    for i in 0..(1 << 23) {
        sqvec.push(i);
    }

    c.bench_function("iter 1", |c| {
        c.iter(|| {
            sqvec.iter().for_each(|x| {
                black_box(x);
            });
        });
    });
}

fn push_benchmark_vec(c: &mut Criterion) {
    // let mut rng = rng();

    c.bench_function("vec push", |c| {
        c.iter(|| {
            let mut vec = Vec::<u32>::new();
            for i in 0..(1 << 23) {
                vec.push(i);
            }
        });
    });
}

criterion_group!(
    benches,
    push_benchmark_sqvec,
    push_benchmark_vec,
    sqvec_iter_bench1
);
criterion_main!(benches);
