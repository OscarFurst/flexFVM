use nalgebra::{Point3, Vector3};

/// Checks if a polygon is convex and its vertices are sorted in a counterclockwise order (right-hand rule).
pub fn check_polygon_convexity(points: &[Point3<f64>]) -> bool {
    let n = points.len();
    if n < 3 {
        return false;
    }
    let mut convex = true;
    let mut normal = Vector3::new(0., 0., 0.);
    for i in 0..n {
        let j = (i + 1) % n;
        let k = (i + 2) % n;
        let v1 = points[i] - points[j];
        let v2 = points[k] - points[j];
        let cross = v1.cross(&v2);
        if i == 0 {
            normal = cross;
        } else if cross.dot(&normal) < 0. {
            convex = false;
            break;
        }
    }
    convex
}

pub fn polyhedron_center(points: &[Point3<f64>]) -> Point3<f64> {
    let mut center = Point3::new(0., 0., 0.);
    for point in points {
        center += point.coords;
    }
    center / points.len() as f64
}

/// Splits a convex polygon into triangles which all share the center of the polygon as a vertex.
pub fn triangulate_polygon(points: &[Point3<f64>]) -> Vec<[Point3<f64>; 3]> {
    let center = polyhedron_center(points);
    let mut triangles = Vec::new();
    for i in 0..points.len() {
        triangles.push([points[i], points[(i + 1) % points.len()], center]);
    }
    triangles
}

pub fn triangle_area(a: &Point3<f64>, b: &Point3<f64>, c: &Point3<f64>) -> f64 {
    let ab = b - a;
    let ac = c - a;
    0.5 * ab.cross(&ac).norm()
}

/// Computes the area of a convex polygon.
pub fn polyngon_area(points: &[Point3<f64>]) -> f64 {
    // start by triangulating the polygon
    let triangles = triangulate_polygon(points);
    triangles
        .iter()
        .map(|[a, b, c]| triangle_area(a, b, c))
        .sum()
}

/// Computes the normal of a convex polygon.
pub fn polygnon_normal(points: &[Point3<f64>]) -> Vector3<f64> {
    let mut normal: Vector3<f64> = Vector3::new(0., 0., 0.);
    for i in 0..points.len() {
        let j = (i + 1) % points.len();
        normal[0] +=
            (points[i].coords.y - points[j].coords.y) * (points[i].coords.z + points[j].coords.z);
        normal[1] +=
            (points[i].coords.z - points[j].coords.z) * (points[i].coords.x + points[j].coords.x);
        normal[2] +=
            (points[i].coords.x - points[j].coords.x) * (points[i].coords.y + points[j].coords.y);
    }
    normal.normalize()
}

/// Computes the volume of a convex polyhedron.
pub fn polyhedron_volume(faces: &[Vec<Point3<f64>>]) -> f64 {
    let mut center = Point3::new(0., 0., 0.);
    for face in faces {
        for point in face {
            center += point.coords;
        }
    }
    center /= (faces.len() * faces[0].len()) as f64;
    let mut volume = 0.;
    for face in faces {
        // triangulate the face
        let triangles = triangulate_polygon(face);
        for [a, b, c] in triangles {
            volume += unsigned_tetrahedron_volume(&center, &a, &b, &c);
        }
    }
    volume
}

fn signed_tetrahedron_volume(
    a: &Point3<f64>,
    b: &Point3<f64>,
    c: &Point3<f64>,
    d: &Point3<f64>,
) -> f64 {
    let ba = a - b;
    let ca = a - c;
    let da = a - d;
    (ba.cross(&ca).dot(&da)) / 6.0
}

fn unsigned_tetrahedron_volume(
    a: &Point3<f64>,
    b: &Point3<f64>,
    c: &Point3<f64>,
    d: &Point3<f64>,
) -> f64 {
    signed_tetrahedron_volume(a, b, c, d).abs()
}
