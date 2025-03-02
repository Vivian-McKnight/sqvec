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
    let mut rng = rng();
    let mut sqvec = SqVec::<u32>::new();
    for _ in 0..(1 << 16) {
        sqvec.push(black_box(rng.random()));
    }

    c.bench_function("iter sqvec", |c| {
        c.iter(|| {
            sqvec.iter_t2().for_each(|x| {
                black_box(x);
            });
        });
    });
}

fn iter_benchmark_vec(c: &mut Criterion) {
    let mut rng = rng();
    let mut vec = Vec::<u32>::new();
    for _ in 0..(1 << 16) {
        vec.push(black_box(rng.random()));
    }

    c.bench_function("iter vec", |c| {
        c.iter(|| {
            vec.iter().for_each(|x| {
                black_box(x);
            });
        });
    });
}

criterion_group!(benches, sqvec_iter_bench1, iter_benchmark_vec);
criterion_main!(benches);
