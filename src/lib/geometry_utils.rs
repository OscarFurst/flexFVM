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

pub fn line_plane_intersection(
    line_point_1: &Point3<f64>,
    line_point_2: &Point3<f64>,
    plane_point: &Point3<f64>,
    plane_normal: &Vector3<f64>,
) -> Option<Point3<f64>> {
    let line_vector = line_point_2 - line_point_1;
    let projection_vector = plane_point - line_point_1;
    let dot_product = plane_normal.dot(&line_vector);
    if dot_product.abs() < f64::EPSILON {
        // CAREFUL: this means the line is parallel to the plane, but it could still be inside the plane!
        return None;
    }
    let t = plane_normal.dot(&projection_vector) / dot_product;
    Some(line_point_1 + t * line_vector)
}

/// Computes the barycentric coordinates of a point in a triangle.
pub fn barycentric_coordinates(
    point: &Point3<f64>,
    a: &Point3<f64>,
    b: &Point3<f64>,
    c: &Point3<f64>,
) -> (f64, f64, f64) {
    let v0 = b - a;
    let v1 = c - a;
    let v2 = point - a;
    let d00 = v0.dot(&v0);
    let d01 = v0.dot(&v1);
    let d11 = v1.dot(&v1);
    let d20 = v2.dot(&v0);
    let d21 = v2.dot(&v1);
    let denom = d00 * d11 - d01 * d01;
    let v = (d11 * d20 - d01 * d21) / denom;
    let w = (d00 * d21 - d01 * d20) / denom;
    let u = 1. - v - w;
    (u, v, w)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_barycentric_coordinates() {
        // values from https://blogs.sas.com/content/iml/2023/06/21/barycentric-coordinates.html
        let a = Point3::new(0.8, -0.2, 0.);
        let b = Point3::new(0.3, 1., 0.);
        let c = Point3::new(0., 0., 0.);
        let point = Point3::new(0.4, 0.6, 0.);
        let (u, v, w) = barycentric_coordinates(&point, &a, &b, &c);
        assert!((u - 0.255814).abs() < 1e-6);
        assert!((v - 0.6511628).abs() < 1e-6);
        assert!((w - 0.0930223).abs() < 1e-6);
    }
}
