//  VoroBox
//
//  Triangulation.swift
//
//  Created by Rob Whitehurst on 4/5/20.
//  Parts of this file are based on Delaunator Copyright (c) 2017, Mapbox
//
// ISC License
//
// Copyright (c) 2017, Mapbox
//
// Permission to use, copy, modify, and/or distribute this software for any purpose
//  with or without fee is hereby granted, provided that the above copyright notice
// and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
// REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
// FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
// INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
// OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
// TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
// THIS SOFTWARE.
//
//  Other parts are
// ISC License
//
//  Copyright © 2020 Rob Whitehurst. oycar@whitehurst.com
// Permission to use, copy, modify, and/or distribute this software for any purpose
//  with or without fee is hereby granted, provided that the above copyright notice
// and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
// REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
// FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
// INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
// OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
// TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
// THIS SOFTWARE.
//
/*:
 ## VoroBox - Conforming Voronoi Diagrams
 
 The underlying triangulation is based on the [Sweep Circle](http://cglab.ca/~biniaz/papers/Sweep%20Circle.pdf) algorithm of
 Biniaz & Dastghaibyfard (2011).
 
 The idea is for points to be included in the triangulation as they enter an expanding circle centred somewhere on the plane.
 ![Sweep Circle Example](./Documentation/Sweep Circle Algorithm.png)
 
 The new triangle is made locally delaunay - through edge flippping, and then once this is done
 the advancing front of the conex hull must be updated. This can be done by searching both
 anti-clockwise and clockwise from the new point, adding extra edges as needed. To facilitate this
 a list of points is maintained sorted by angle from an initial point inside  the convex hull.
 
 This initial seed point need not be the same as the radial origin, but can be if the radial point is
 within the convex hull. To [ensure this](https://cs.stackexchange.com/questions/63362/delaunay-sweep-circle-initialization)
 a good seed point is one is to take a point from the set, its nearest neighbour, and a 3rd point that creates the smallest circumcircle with them will always work
 as this is a delaunay triangle (nearest-neighbours are always delaunay edges, and the 2 triangles based on this edge have to create an empty circumcircle, so
 one of them has to be the smallest) which makes the circumcircle empty of other points, so using its origin as the radial origin is correct.
 But (in principle) the circumcentre might be outside the convex hull; so using the incircle centre would work to for the radial origin & the angular origin.
 
 
 */



//
import Foundation
import Files


// Define statics
internal let Epsilon:Double = pow(2.0, -53)
// Shewchuk error bounds
let squaredThreshold:Double = (2.0 +        Epsilon) * Epsilon
internal let circleThreshold:Double = (10.0 + 96.0 * Epsilon) * Epsilon
internal let orientThreshold:Double =  (3.0 + 16.0 * Epsilon) * Epsilon

internal let ProjectFolder = ProcessInfo.processInfo.environment["HOME"]!  + "/CloudStation/Projects/VoroBox"
internal let ZoneFolder = ProjectFolder + "/Zones"
internal let OutputFolder = ProjectFolder + "/Output"
internal var showCount = 0
internal let showLimit = 40

// Codes with magnitudes less than BoundaryFlag are available for special purposes
internal let BoundaryEdge = -1 // On the expanding hull
internal let EmptyEdge = -2 // Useful for removed edges
internal let UnusedEdge = -3
internal let NoEdge = -4
internal let FinishedEdge = -5
internal let EmptyVertex = -888

// Vertex ordering - HSPQR
// Last vertex R is unique to any corner and always present in the edge loops
internal let hIndex = 0, sIndex = 1, pIndex = 2, qIndex = 3, rIndex = 4

// Bounding Box
internal var boundingBox:Array<Array<Double>> = [[Double.infinity, Double.infinity], [-Double.infinity, -Double.infinity]]

 enum triangulationError: Error {
    case initError(String)
 }

// The main structure
struct Triangulation: Codable {
  // Store all trangulation coordinates
  // The count of points, their co-ordinates
  static var coords = Array<Double>()
  static var properties = Array<Int>() //  Index of properties associated with this vertex
  static var pointCount = 0
  
  // The singleton triangulation
  static var triangulation = Triangulation()
  
  var vertices = Array<Int>()
  var halfEdges = Array<Int>()
  var emptyTriangles = Set<Int>()

  var numberEdges: Int      = 0
  
  // The edge stack
  static var edgeStack = Array<Int>()
    
  // Hull properties (too many?)
  // The triangulation boundary
  static var hullNext = Dictionary<Int, Int>()
  static var hullPrev = Dictionary<Int, Int>()
  
  // These are loops of vertices
  static var loopVertices = Dictionary<Int, Array<Int>>()
  static var loopHooks = Dictionary<Int, Int>() // loopHook[p] => [q]; edge defined by p->q

  // Convex hull properties
  var convexHullEdges = Dictionary<Int, Int>()
    
  // A list of Voronoi zones and the associated vertex index
  var zoneList = Array<Array<Int>>()
  var zoneVertex = Array<Int>()
  
  // Need information to construct voronoi cells
  var triangleCentres = Dictionary<Int, Array<Double>>()
  var triangleCount = 0
}

extension Triangulation {
  init(size numberPoints:Int) {
    // arrays that will store the triangulation
    let size = 3 * max(2 * numberPoints - 5, 0)
    vertices = Array(repeating:0, count:size)
    halfEdges = Array(repeating:UnusedEdge, count:size)
  }
  
  mutating func growTriangulation(to size:Int) {
    let existingSize:Int = vertices.count
    let extraSize = 3 * max(2 * size - 5, 0) - existingSize
    if extraSize > 0 {
      // Grow the two arrays
      vertices.append(contentsOf:Array(repeating:0, count: extraSize))
      halfEdges.append(contentsOf:Array(repeating:UnusedEdge, count: extraSize))
    }
  }
  
  // Add a set of vertices
  mutating func addVertices(indices zoneIndices: Array<Int>) throws ->  Dictionary<Int, Int> {
    var hashSize = 0
    var hashFactor: Double = 0

    // Get a hash key - slightly simplified
    func hashKey(_ dx:Double, _ dy:Double) -> Int {
      // Pseudo angle phi
      var phi = dx / (abs(dx) + abs(dy))
      phi = dy > 0 ? 3 - phi : 1 + phi // Ranges from [0..4]
      return Int(hashFactor * phi.rounded(.down)) % hashSize
    }
    
    // Use the number of local indices...?
    let numberPoints = zoneIndices.count
    var convexHullPrev: Dictionary<Int, Int> = [:] // edge to prev edge
    var convexHullNext: Dictionary<Int, Int> = [:] // edge to next edge
    
    // Distances from the origin plus index ids
    var dists = Array(repeating:0.0, count:numberPoints)
    var ids   = Array(repeating:0, count:numberPoints)
    
    // Get a random vertex
    let i0 = zoneIndices.randomElement()!
    let i0x = Triangulation.coords[2 * i0]
    let i0y = Triangulation.coords[2 * i0 + 1]
    
    // Set the ids array and the distances from this vertex
    for (i, j) in zoneIndices.enumerated() {
      // Re-index the points
      // ids are local ids
      ids[i] = i
      if j != i0 {
        dists[i] = dist(Triangulation.coords[2 * j], Triangulation.coords[2 * j + 1], i0x, i0y)
      } else { dists[i] = 0 }
    }
    
    // sort the ids by distance from this vertex
    ids.sort { dists[$0] < dists[$1] }
    
    // find the point closest to the seed
    var i1 = zoneIndices[ids[1]]
    var i1x = Triangulation.coords[2 * i1]
    var  i1y = Triangulation.coords[2 * i1 + 1]
    
    // We need a seed triangle
    // So compute the circumcircle radius of each point
    // and save that with the minimum radius
    
    // We now construct a candidate seed triangle
    var r2min = Double.infinity
    var i2 = 2
    for k in 2..<zoneIndices.count {
      let j = ids[k]
      // Once this point is more than double the minimum radius
      // away from point i0 the circumradius cannot get any smaller
      if dists[j] > 4.0 * r2min { break }
      
      //
      let j2 = zoneIndices[j]
      let j2x = Triangulation.coords[2 * j2]
      let j2y = Triangulation.coords[2 * j2 + 1]

      // The circumradius
      let r2 = circumRadius(i0x, i0y, i1x, i1y, j2x, j2y)
      if r2 < r2min {
        r2min = r2
        
        // Save this index
        i2 = j2
      }
    }
    
    var i2x = Triangulation.coords[2 * i2]
    var i2y = Triangulation.coords[2 * i2 + 1]
    
    // Co-linear?
    if r2min == Double.infinity {
      throw triangulationError.initError("All vertices are co-linear, no triangulation possible")
    }
    
    // swap the order of the seed points
    // to ensure seed triangle has
    // an anticlockwise orientation
    if orient(i0x, i0y, i1x, i1y, i2x, i2y) < 0 {
      // Swap i1 & i2
      (i1, i2) = (i2, i1)
      (i1x, i2x) = (i2x, i1x)
      (i1y, i2y) = (i2y, i1y)
    }
    
    //
    // We have a seed triangle
        
    // The initial triangle defines the hull
    convexHullEdges[i0] = numberEdges
    convexHullEdges[i1] = numberEdges + 1
    convexHullEdges[i2] = numberEdges + 2
    
    // convexHullNext[i] -> points to next index on the hull in an anti-clockwise direction
    // convexHullPrev[i] -> points to next index on the hull in an      clockwise direction
    convexHullNext[i0] = i1; convexHullPrev[i2] = i1
    convexHullNext[i1] = i2; convexHullPrev[i0] = i2
    convexHullNext[i2] = i0; convexHullPrev[i1] = i0
    
    // Add the initial triangle
    let _ = addTriangle(i0, i1, i2, BoundaryEdge, BoundaryEdge, BoundaryEdge)

    // Most work done when no internal points
    if numberPoints > 0 {
      // The size allowed for the hash array ordered by angle
      // Since it is for  the hull it is approx. the square root of the total
      // number of points - but not if all the points are on the perimeter..
      hashFactor = Double(numberPoints).squareRoot().rounded(.up)
      hashSize = Int(hashFactor)
      hashFactor *= 0.25
      var hullHash = Array(repeating:-1, count:hashSize) // angular edge hash
      
      // Get the circumcentre of the seed triangle
      // No points should be inside the
      let seed:Array<Double> = circumCentre(i0x, i0y, i1x, i1y, i2x, i2y)
      
      // Set distances
      // The indices are sparse - so pack them into tighter arrays
      for (i, j) in zoneIndices.enumerated() {
        ids[i] = i
        dists[i] = dist(Triangulation.coords[2 * j], Triangulation.coords[2 * j + 1], seed[0], seed[1])
      }

      // sort the points by distance from the seed triangle circumCentre
      ids.sort {dists[$0] <  dists[$1]}
      
      // The hash is based on a pseudo-angle so it increases with azimuth
      // Should be able to use a dictionary
      // Sort keys in angle order
      // Stop when next key is desired one (as in mpx sorting)
      hullHash[hashKey(i0x - seed[0], i0y - seed[1])] = i0
      hullHash[hashKey(i1x - seed[0], i1y - seed[1])] = i1
      hullHash[hashKey(i2x - seed[0], i2y - seed[1])] = i2
      
      // A list to hold boundary flips
      var boundaryFlips = Array<Array<Int>>?([])

      // Now start building the triangulation
      // The previous point
      var xp:Double = 0, yp:Double = 0
      
      // A point to start looking
      var start = i0

      // Labels! Just like good old FORTRAN IV
      projectPoints: for (_, index) in ids.enumerated() {
        func patchHull(with flips:Array<Array<Int>>) {
          // Patch the hull
          for f in flips {
            // Patch this edge
            convexHullEdges[vertices[f.first!]] = f.first!
          }
        }
        
        let i = zoneIndices[index]
        let x = Triangulation.coords[2 * i]
        let y = Triangulation.coords[2 * i + 1]
        
        // skip near-duplicate points
        if isZero(x - xp) && isZero(y - yp) {
          continue projectPoints
        }

        // Save previous point
        xp = x
        yp = y
        
        // skip seed triangle points
        if (i == i0 || i == i1 || i == i2) {continue projectPoints}
      
        // find a visible edge on the convex hull using edge hash
        // This is the projection step in the algorithm
        // Get the current pseudo-angle hash key
        let key = hashKey(x - seed[0], y - seed[1])
        
        // The projection of the point on the convex hull
        // will cross an edge - find the start of this edge
        // This can have hash collisions and so miss some
        // edges - also stale entries x are present
        // which are labelled by setting convexHullNext[x] => x
        for j in 0..<hashSize {
          // Loop over all the hash keys
          start = hullHash[(key + j) % hashSize]
          
          // Eventually get to a non-empty hash entry
          // Then check the next entry - if it is identical keep going because
          // that means the point was removed
          // When one is found it is the next point on the hull with a greater azimuth
          if (start != -1 && start != convexHullNext[start]) {break}
        }
        
        // Get to here
        // get the previous point - we have bracketed  the projected point
        // azimuth(e) <= azimuth(i) <= azimuth(q)
        // - but not sometimes due to collisions
        var e = convexHullPrev[start]!
        var q = start
        start = e
        
        // Check for positive (anti-clockwise) orientation of the triangle [i, q, e] where the projected point is i
        // So if [i, q, e] is clockwise this test fails and the triangle is replaced in the loop
        // This can happen with hash collisions
        // Another problem with the hash is the
        while (orient(x, y,
                       Triangulation.coords[2 * q], Triangulation.coords[2 * q + 1],
                       Triangulation.coords[2 * e], Triangulation.coords[2 * e + 1]) <= 0) {
                        if q == start {
                          // Can't add this point! INSIDE convex hull
                          throw triangulationError.initError("Found vertex \(i) [\(Triangulation.coords[2*i]), \(Triangulation.coords[2*i+1])] inside convex hull")
                        }
                        
                        // Move around  the hull by one triangle until [i, q, e] is anti-clockwise
                        //
                        //       i            i                 i
                        //       |           /|                /|
                        //       |          / |               / |
                        // n <-- q     =>  n--q  relabel =>  q--e
                        //       |
                        //       |
                        //       e <--
                        //   +ve <== Direction
                        e = q
                        q = convexHullNext[e]!
        }
        
        
        // add the first triangle from the point    i
        // This is actually just [e, i, q]         / \
        //                                        q---e
        //
        // The linked edges are BoundaryEdge (outside the hull) for (e->i) and (i->q) and hillTri[e] for (q->e)
        var t = addTriangle(e, i, q, BoundaryEdge, BoundaryEdge, convexHullEdges[e]!)
        
        // Number of edges has increased by 3; t is the old number of vertices
        // Make the new triangle locally delaunay
        convexHullEdges[i] = try! flipEdges(t + 2, list: &boundaryFlips)
        convexHullEdges[e] = t  // keep track of boundary vertices on the hull

        // Patch the hull
        patchHull(with: boundaryFlips!)
        boundaryFlips!.removeAll(keepingCapacity: true)
        
        // We added one triangle - so the hull size goes up by one
        //
        //     i
        //    / \
        //   q---e        The hull edge (e->q) has been replaced by two new edges (e->i) and (i->q)
        //    \ /
        //     p
        //
        // Replaced convexHullEdges[e] and added convexHullEdges[i]
        // Notice that convexHullEdges has lots of zero entries - its size is not being changed here, just some zeros are being replaced
        // It might be better if they were initialized with a negative number
        
        // Move forward around the hull
        //
        //        i             i
        //       / \           / \
        //      /   \         /   \
        //    q ---- n   <=  q ---- e
        //        Anticlockwise
        var n = q
        q = convexHullNext[q]!
        
        // Walk around anti-clockwise; add the triangle [n, i, q]
        //
        //
        // This orientation test is on the triangle [n, i, q] => true if anticlockwise
        while (orient(x, y,
                      Triangulation.coords[2 * q], Triangulation.coords[2 * q + 1],
                      Triangulation.coords[2 * n], Triangulation.coords[2 * n + 1]) > 0) {
                        
                        // Add this triangle to the triangulation
                        t = addTriangle(n, i, q, convexHullEdges[i]!, BoundaryEdge, convexHullEdges[n]!)
                        
                        // We added edge i=>q
                        // Was it already there from an earlier zone?
                        
                        // Fix locally non-delaunay edges
                        // Returns hull edge i on return (the new edge)
                        convexHullEdges[i] = try! flipEdges(t + 2, list: &boundaryFlips)
                        
                        // Patch the hull
                        patchHull(with: boundaryFlips!)
                        boundaryFlips!.removeAll(keepingCapacity: true)
                        
                        // This indicates the edge n is no longer on the hull
                        convexHullNext[n] = n // mark as removed
                        
                        // Continue going anti-clockwise
                        n = q
                        q = convexHullNext[n]!
        }
        
        // walk backward from the other side, adding more vertices
        // This is a clockwise search
        // Only carried out when the projected point i split
        // the edge e->q to create an anti-clockwise triangle [e, i, q]
        // Because otherwise we have
        //
        //       i
        //       |
        //       |
        // n <-- q
        //       |
        //       |
        //       e <--
        //   +ve <== Direction
        //
        // And so the point e is not visible from i
        if (e == start) {
          
          // Move backwards around the hull
          //
          //        i             i
          //       / \           / \
          //      /   \         /   \
          //    q ---- e   =>  e ---- q
          //         Clockwise
          q = convexHullPrev[e]!
          
          // This orientation test is on the triangle [i, e, q] => true if anticlockwise
          while (orient(x, y,
                        Triangulation.coords[2 * e], Triangulation.coords[2 * e + 1],
                        Triangulation.coords[2 * q], Triangulation.coords[2 * q + 1]) > 0) {
                          
                          // Add this triangle to the triangulation
                          t = addTriangle(q, i, e, BoundaryEdge, convexHullEdges[e]!, convexHullEdges[q]!)
                          
                          // We added edge q=>i                          
                          // Fix locally non-delaunay edges
                          // Ignore return value (which is edge i - not on the hull)
                          _ = try! flipEdges(t + 2, list: &boundaryFlips)
        
                          // Patch the hull
                          patchHull(with: boundaryFlips!)
                          boundaryFlips!.removeAll(keepingCapacity: true)
                      
                          // Add new hull edge q
                          convexHullEdges[q] = t
                          
                          // This indicates the edge e is no longer on the hull
                          convexHullNext[e] = e // mark as removed
                          
                          // Continue going clockwise
                          e = q
                          q = convexHullPrev[e]!
          }
        }
        
        // update the hull
        // These are vertices
        convexHullPrev[i] = e ; convexHullNext[e] = i
        convexHullPrev[n] = i ; convexHullNext[i] = n
        
        // save the two new edges in the hash table
        hullHash[hashKey(x - seed[0], y - seed[1])] = i
        hullHash[hashKey(Triangulation.coords[2 * e] - seed[0], Triangulation.coords[2 * e + 1] - seed[1])] = e
      }
    }
    
    // Return the convex hull of vertices
    return convexHullNext
  }

  private mutating func link(_ a: Int, _ b: Int) {
    // Just link half edges
    halfEdges[a] = b
    if b > BoundaryEdge {
      halfEdges[b] = a
    }
  }

  // add a new triangle given vertex indices and adjacent half-edge ids
  private mutating func addTriangle(_ i0: Int, _ i1: Int, _ i2: Int,
                            _  a: Int, _  b: Int, _  c: Int) -> Int {

    // Save the value of numberEdges
    let t = numberEdges

    // Add the triangle [i0, i1, i2]
    // These actually indicate halfEdge ids
    vertices[t] = i0
    vertices[t + 1] = i1
    vertices[t + 2] = i2

    // Connect the paired halfEdges
    link(t, a)
    link(t + 1, b)
    link(t + 2, c)

    // Added three extra halfEdges
    numberEdges += 3
    
    // Return the original number of halfEdges
    return t
  }

  // Hull functions
  mutating func join(loop convexHullNext:Dictionary<Int, Int>) throws {
    // Join convexZone hulls together
    // Strip empty values from the hull
    var convexNext = Dictionary<Int, Int>()
    var convexPrev = Dictionary<Int, Int>()
    for (p, q) in convexHullNext {
      if p != q {
        convexNext[p] = q
        convexPrev[q] = p
      }
    }
    
    // Not empty
    if convexNext.isEmpty { return }
    //if showMe > 1  {
      //print("Join Hulls")
      //showVoronoi(label: "jl_", type: "Delaunay")
    //}
    
    // Save any matched edges
    var matchedEdges = Dictionary<Int, Int>()
    
    // Get each edge from the convex hull
    let hullStart = convexNext.first!.value
    var p:Int
    
    // Most efficient to search the existing hull loops before adding any new edges
    if !Triangulation.hullNext.isEmpty {
      //         We have the vertices p, q
      //         And the edges        e, f
      //         p -> e -> q -> f
      //
      p = hullStart
      var e = convexHullEdges[p]!
      
      // Get each convex hull edge
      convexLoop: repeat {
        // Next vertex
        let q = convexNext[p]!
        
        // Before storing this in the triangulation loops
        // check that it is not a matching edge
        // Since each edge is unique we
        // look for matching vertices
        //
        //   \ d       f /
        //    \    e    /
        //     p-------q
        //    /    b    \
        //   / c       a \
        //
 
        hullLoop: for (b, c) in Triangulation.hullNext {
          // Only gets to  here if hullNext is not empty
          // The convex hull   is edge e p->q
          // The triangulation is edge b q->p
          if p == vertices[c] && q == vertices[b] {
            // This is a match
            // Order does matter - matchedEdges[Convex] => Triangulation
            matchedEdges[e] = b
            break hullLoop
          }
        }
        
        // Get next convex hull edge
        e = convexHullEdges[q]!
        
        // Next vertex
        p = q
      } while hullStart != p
    }
  
    // Now add the convex edges to the loops
    // Don't join the zones yet because if
    // the zones are flipped the loop edges can change
    p = hullStart
    var e = convexHullEdges[p]!
    var q:Int, f:Int
    repeat {
      q = convexNext[p]!
      f = convexHullEdges.removeValue(forKey: q)!

      // Connect loops
      Triangulation.hullNext[e] = f
      Triangulation.hullPrev[f] = e
      
      // Move around
      p = q
      e = f
    } while hullStart != p
    
    // Finally join the matched edges
    // Edges can change!
    // So peel off the first edge one at a time
    while !matchedEdges.isEmpty {
      let firstMatch = matchedEdges.first!
      let e = firstMatch.key
      let b = matchedEdges.removeValue(forKey: e)!
      
      // The vertices 
      p = vertices[e]
      q = vertices[b]
      
      //
      //      convex
      //   \ d       f /
      //    \    e    /
      //     p-------q
      //    /    b    \
      //   / c       a \
      //   Triangulation
      //
      
      let a = Triangulation.hullPrev[b]!, f = Triangulation.hullNext[e]!
      let c = Triangulation.hullNext[b]!, d = Triangulation.hullPrev[e]!
      
      // New entries
      Triangulation.hullPrev[f] = a; Triangulation.hullNext[a] = f
      Triangulation.hullPrev[c] = d; Triangulation.hullNext[d] = c
      
      // Clear old pointers
      Triangulation.hullPrev[e] = nil; Triangulation.hullNext[e] = nil
      Triangulation.hullPrev[b] = nil; Triangulation.hullNext[b] = nil

      // Join the halfEdges
      halfEdges[e] = b
      halfEdges[b] = e

      // four-way flip -> should become two way flips
      //
      // This flips into the triangulation first
      var flipBoundaries = Array<Array<Int>>?([])
      _ = try! flipEdges(b, list: &flipBoundaries, flipFour: true)
      
      // Patch the loops
      // In order!
      for f in flipBoundaries! {
        // Splice new edge y for old edge x
        //
        //         x
        //     p-------q
        //    /    y    \
        //   / c       a \
        //   Triangulation
        //
        let new = f.first!, old = f.last!

        // Sometimes an edge can be swapped out
        loopEdges(edge: new, replaces: old)

        // Patch the matchedEdges
        // Is this correct?
        if let value = matchedEdges[old] {
          matchedEdges[new] = value
          matchedEdges[old] = nil
        } else if let element = matchedEdges.first(where: {$0.value == old }) {
          matchedEdges[element.key] = new
        }
      } // End of each flip boundary
    } // End of matched
    // Finished
  }
  
  // Need to identify the hull loops
  mutating func hullLoops() -> Int {
    // We just want a list of all the vertices that are on the hull
    // The only complication is that there is more than one loop
    // The linked list of edges
    var hullEdges = Triangulation.hullNext
    
    var index = 0
    var hookIndex = 0
    while !hullEdges.isEmpty {
      // When reproducible is set get the vertex with ab extremum value of x
      var loopStart:Int
      
      // Only do this when debugging
      if showMe > 0 {
        let x = hullEdges.sorted {
          let a = vertices[$0.key]
          let b = vertices[$1.key]
          let ax = Triangulation.coords[2 * a]
          let bx = Triangulation.coords[2 * b]
          if ax == bx {
            let ay = Triangulation.coords[2 * a + 1]
            let by = Triangulation.coords[2 * b + 1]
            return ay >= by
          }
        
          // Get the relative orientation
          return ax > bx
        }
        let y = x.first!.key
        
        loopStart = y
      } else {
        loopStart = hullEdges.first!.value
      }
      
      // First edge
      var e:Int? = loopStart
      repeat {
        // Save this vertex
        //Triangulation.loopVertices[index] = [vertices[e!]]
        index += 1
        
        // Next edge
        e = hullEdges.removeValue(forKey: e!)
      } while loopStart != e
      
      // One loop complete
      Triangulation.loopHooks[hookIndex] = loopStart
      hookIndex = index
    }
    
    // All done
    return Triangulation.loopHooks.count
  }
  
  // Conforming Voronoi
  // First the hubs need to defined
  mutating func voronoiConforming() throws {

    // First pass finds the generating points to are used to obtain the constructors
    var hubPoints = Array<Array<Double>>()
    var hubRadii  = Dictionary<Int, Double>()
    var isLinear  = Set<Int>()
    if showMe != 0 {
      print("Pass 2")
      print("\tFind hub generators")
      showVoronoi(label: "fh_", type: "Delaunay")
    }
    
    // Get each hull loop
    let loopKeys = Triangulation.loopHooks.keys.sorted()
    for iº in loopKeys {
      // Get a boundary edge from this loop
      let loopStart = Triangulation.loopHooks[iº]!
      
      // Find the first vertex
      var eº = Triangulation.hullPrev[loopStart]!
      var fº = loopStart
      var h  = vertices[eº]
      var j  = vertices[loopStart]
      var hx = Triangulation.coords[2 * h], hy = Triangulation.coords[2 * h + 1]
      var jx = Triangulation.coords[2 * j], jy = Triangulation.coords[2 * j + 1]
      
      // Set Size
      var hSize:Double = 0
      var firstSize:Double? // Initially zero
      
      // Get the next index and hub
      roundLoop: repeat {
        // Update Edges
        eº = fº
        fº = Triangulation.hullNext[fº]!
        
        // Update vertices
        let g = h
        h  = j
        j = vertices[fº]

        // Save the previous size 
        let gSize = hSize
        
        // The points associated with these vertices
        let gx = hx
        let gy = hy
        hx = jx
        hy = jy
        jx = Triangulation.coords[2 * j]
        jy = Triangulation.coords[2 * j + 1]
        
        // We need a vertex plus the two edges
        //
        //                  g
        //                 /
        //                /
        //               /
        //              h ::::::::: <angle bisector> ::::
        //               \
        //                \
        //                 \
        //                  j
        //
        // The triangle incentre lies on the angle bisector
        
        // The sense of the triangle [ghj]
        // If sense is positive a convex corner
        // If sense is zero not a corner - colinear
        var sense = orient(gx, gy, hx, hy, jx, jy)
        
        // Now we need incircle Centre or (for colinear points) this is just the hub
        // All colinear points will have at least one true triangle
        // which is adjacent to it so a reliable result for the
        // location of the generator point can be derived
        var centre:Array<Double>? = nil
        if 0 != sense {
          centre = inCentre(gx, gy, hx, hy, jx, jy)
          
          // Provide a limiting value for the hub size
          let dx = centre![0] - hx
          let dy = centre![1] - hy
          let r2 = dx * dx + dy * dy
          
          // Check this
          if r2 < squaredThreshold {
            // Treat as colinear
            // First though treat the centre as if it were convex
            if sense < 0 {
              centre = [hx - dx, hy - dy]
            }
            sense = 0
          }
        }
        
        // The hSize
        hSize = hubRadii[h, default:Double.infinity]
        
        // Is this colinear?
        if 0 == sense { isLinear.insert(hubPoints.count) }
        
        // Now loop over the spokes connected to the hub h
        // This is an anticlockwise loop
        // Ends on the outgoing edge from vertex (g)
        
        var a = eº, a2:Int
        repeat {
          // Get this triangle
          let a0 = 3 * (a / 3)
          let a1 = a0 + (a + 1) % 3
          a2 = a0 + (a + 2) % 3
          
          // Update internal vertices (we can be more efficient than this if it works)
          let i1 = vertices[a1]
          let i2 = vertices[a2]
          
          // Coordinates
          let i1x = Triangulation.coords[2 * i1]
          let i1y = Triangulation.coords[2 * i1 + 1]
          let i2x = Triangulation.coords[2 * i2]
          let i2y = Triangulation.coords[2 * i2 + 1]
          
          // if i1 is j & i2 is g this is done (only one triangle ie [hjg]
          // Only one triangle never occurs for a colinear point
          let c = j == i1 && g == i2 ? centre! : inCentre(hx, hy, i1x, i1y, i2x, i2y)
          
          // Get r2
          let dx = c[0] - hx
          let dy = c[1] - hy
          let r2 = dx * dx + dy * dy
          
          // We need the incircle radius for this triangle
          if r2 < hSize {
            // Note this
            hSize = r2
          }
          
          // Next edge
          a = halfEdges[a2]
          
          // Finished if a boundary
        } while BoundaryEdge < a

        // For reflex corners the size of the preceeding hub is important 
        var hubSize = hSize
        if eº != loopStart {
          if sense <= 0 {
            // Treat co-linear as reflex? 
            if gSize < hSize {
              hubSize = gSize
            } else if hSize < hubRadii[g]! {
              // Correct previous entry 
              hubRadii[g] = hSize
            }
          }
        } else if sense <= 0 { 
          // Save first 
          firstSize = hSize
        }

        // Save the hSize and the incentre at this index
        if sense == 0 {
          // We need to obtain a point perpendicular
          // to the line connecting g & j at the hub h of length squareRoot(hubSize)
          let r = perpendicular(scale: hubSize, hx, hy, gx, gy, jx, jy)
          hubPoints.append(r)
        } else {
          hubPoints.append(centre!)
        }
     
        // The hub radii is accumulated for every hub h -
        // so check this hub was not recorded before on another loop
        if let r2 = hubRadii[h] {
          if hubSize < r2 {
            // Replace the entry
            hubRadii[h] = hubSize
          }
        } else { // New entry
          hubRadii[h] = hubSize
        }
        
        // Debugging
        if showMe > 1 {
          print("\thub => \(h)")
          print("\tindex => \(hubPoints.count - 1)")
          print("\tcentre => \(hubPoints.last!)")
          print("\tradius => \(hubRadii[h]!.squareRoot())")
        }
        
      } while loopStart != fº

      // Still have not checked the first size 
      // Here last vertex on loop is h, first is j 
      if let s = firstSize {
        // Only non nil when not convex
        let g = h
        h = j 

        let gSize = hSize
        hSize = s 

        // Both might need setting 
        if gSize < hubRadii[h]! {
          hubRadii[h] = gSize
        } 
        
        if hSize < hubRadii[g]! {
          // Correct previous entry 
          hubRadii[g] = hSize
        }
      } 


    } // End of each loop

    // Third pass computes the constructor vertices
    // Get the vertices used to construct the conforming Voronoi diagram
    // var hubIndex = 0
    if showMe != 0 {
      print("Pass 3")
      print("\tCompute hub sentinels")
    }
    
    // Copy the loop hooks
    // but later ones might be scrambles by removeVertex
    var loopHooks = Triangulation.loopHooks
    
    // Delete the old copy - because new keys might clobber old ones
    Triangulation.loopHooks.removeAll(keepingCapacity: true)
    
    // Get each hull loop
    // for jº in loopHooks.keys.sorted() {
    for jº in loopKeys {
      // Clean loop hooks
      var iº = jº
      let loopStart = loopHooks[jº]!
      
      var eº = Triangulation.hullPrev[loopStart]!
      var fº = loopStart
      var h  = vertices[eº]
      var j  = vertices[loopStart]
      var hx = Triangulation.coords[2 * h], hy = Triangulation.coords[2 * h + 1]
      var jx = Triangulation.coords[2 * j], jy = Triangulation.coords[2 * j + 1]
            
      // Get the next index and hub
      roundLoop: repeat {
        //
        // The next step is to generate the required new points
        // There are four possible sentinel points
        // But one is always redundant, and for colinear hubs two are redundant
        //
        var pº = EmptyVertex, qº = EmptyVertex, rº = EmptyVertex, sº = EmptyVertex
        
        // Update Edges
        eº = fº
        fº = Triangulation.hullNext[fº]!
        
        // Update vertices
        h  = j
        j = vertices[fº]
        
        // The points associated with these vertices
        let gx = hx
        let gy = hy
        hx = jx
        hy = jy
        jx = Triangulation.coords[2 * j]
        jy = Triangulation.coords[2 * j + 1]
        
        // Now we need the incircle at [ghj]
        // In some cases this is not actually an incircle centre
        let centre = hubPoints[iº]
        
        // Compute the hub centre
        var dx = centre[0] - hx
        var dy = centre[1] - hy
        let r2 = dx * dx + dy * dy
        
        // The hub radius squared
        let h2 = hubRadii[h]!
        
        // Scale the distance from hx to the centre by the
        // size of the hub radius
        let rd = 0.5 * (h2 / r2).squareRoot()
        //let rd = 0.25 * (h2 / r2).squareRoot()

        // Now finally we can compute the location
        dx = hx + rd * dx
        dy = hy + rd * dy
        
        // The generator point is noted
        Triangulation.coords.append(dx)
        Triangulation.coords.append(dy)

        // Sentinels pick up generator codes 
        Triangulation.properties.append(Triangulation.properties[h])
                
        // Increment the point count
        let d = Triangulation.pointCount
        Triangulation.pointCount += 1
        
        // We have a generator point d on the bisector; internal for a convex hub; external for a reflex hub
        //
        //
        //                Convex Hub                      Reflex Hub
        //
        //                      g                            g
        //                   . /                             .\
        //            sº    . /                               .\   pº
        //                 . /                                 .\
        //       ---------- h ---pº-----            ------rº---- h --------
        //                 . \                                 ./
        //            rº    . \                               ./   qº
        //                   . \                             ./
        //                      j                            j
        //
        //      External Vertices rº & sº          Internal Vertices pº & qº
        //
        //
        //  (For a colinear point pº == qº AND rº == sº
        
        
        
        if !isLinear.contains(iº) { // hubIndex
          // This is an internal vertex for a convex corner
          // But an external vertex for a reflex corner
          let sense = orient(gx, gy, hx, hy, jx, jy)
          
          // The reflected points
          let d_out = reflect(dx, dy, hx, hy, jx, jy)
          let d_in  = reflect(dx, dy, gx, gy, hx, hy)
          
          // The properties of h are copied - twice
          Triangulation.properties.append(Triangulation.properties[h])
          Triangulation.properties.append(Triangulation.properties[h])
          
          if sense > 0 {
            // Internal generator on convex corner
            pº = d
            qº = pº
            
            // Convex corner
            // Two External points
            // Inward reflexion
            sº = d_in
            
            // Outward reflexion
            rº = d_out
          } else {
            // External generator on reflex corner
            rº = d
            sº = rº
            
            // Reflex corner
            // Need the *internal* points here
            // Inward reflexion
            pº = d_in
            
            // Outward reflexion
            qº = d_out
          }
        } else {
          // Colinear point
          
          // The properties of h are copied - once
          Triangulation.properties.append(Triangulation.properties[h])
          
          // Internal generator same as on a convex corner
          pº = d
          qº = pº
          
          // External generator same as on a reflex corner
          sº = reflect(dx, dy, gx, gy, jx, jy)
          rº = sº
        }
        
        // Next hub index
        // hubIndex += 1
        
        // Save the vertices
        // First is the hub; last is always either the hub or rº
        Triangulation.loopVertices[iº] = [h, sº, pº, qº, rº]
        
        if showMe > 1 {
          print("\thub => \(h)")
          print("\t\tindex => \(iº)")
          print("\t\tlist => \(Triangulation.loopVertices[iº]!)")
        }
        
        
        // Save first encounter as the loop hook - need to maintain the
        // ordering of the loops; currently this is g->h
        if jº == iº {
          // The first key will refer to s
          Triangulation.loopHooks[sº] = j
        }
        
        // Next index
        iº += 1
      } while loopStart != fº
    } // End of Third Pass
    
    // Third pass folds out loops; removing the hub vertices
    let numberHubs:Double = Double(Triangulation.loopVertices.count)
    let hubStep:Int = max(Triangulation.loopVertices.count / 100, 10)
    var countHubs = 0, hubThreshold = hubStep
    if showMe > 0 {
      print("Pass 4")
      print("\tFold out hubs")
      showVoronoi(label: "fo_", type: "Delaunay")
    }
    
    // This adds the hub vertices and adjusts the triangulation
    // vertices are not changed but edges and hooks to edges are
    // Grow the triangulation
    growTriangulation(to: Triangulation.pointCount)
    
    // We need to reindex the loops' vertices
    let loopVertices = Triangulation.loopVertices
    Triangulation.loopVertices.removeAll(keepingCapacity: true)
    
    // Get each hull loop
    // Need to update the loopHooks again!
    loopHooks = Triangulation.loopHooks
    Triangulation.loopHooks.removeAll(keepingCapacity: true)
    
    // Each loop hook must be obtained in the correct order...
    for (key, j) in loopHooks {
      // The key is unique to a given list of vertices
      // find it
      var iº = loopVertices.first(where: {key == $0.value[sIndex]})!.key
      var list = loopVertices[iº]!
      
      // Colinear hubs - all the generator points are collapsed to h
      // Find the edge joining the vertices s & r
      let loopStart = findLoopEdge(from: list[hIndex], to: j)
      
      // Set the initial edge to start at h
      // Edge f goes from h->j
      var f = loopStart
      var e = Triangulation.hullPrev[f]!
      
      // So the stop vertex is g or 'loopStop'
      var loopStop = EmptyVertex // Not initialized yet
            
      // Get the next index and hub
      var hookVertex = EmptyVertex
      foldOut: repeat {
        // Get the vertices waiting to be added
        iº += 1

        // Unpack the vertices
        let pº = list[pIndex], qº = list[qIndex]
        let rº = list[rIndex], sº = list[sIndex]

        // Reindex these vertices
        Triangulation.loopVertices[rº] = list
        
        // Adjust the stop
        if EmptyVertex == loopStop {
          loopStop = sº
        }
        
        // Update loop hook
        if EmptyVertex == hookVertex {
          // This is the new hook for this loop
          hookVertex = rº
          Triangulation.loopHooks[sº] = rº
        }
        
        // Debugging 
        let foldLabel = "fo_" // "fo_\(list[hIndex])_"

        // Convex hubs always look like this (labels refer to internal edges),
        //    eº and fº point in the opposite direction to e & f
        //
        //                Convex Hub
        //            eº+1
        //         sº----- g
        //          \     /
        //      eº+2 \ eº/e
        //            \ /
        //             h ----- pº
        //            / \
        //      fº+1 / fº\f
        //          /     \
        //         rº ---- j
        //            fº+2
        //
        
        // Is this a convex hub?
        showVoronoi(label: "fo_", type: "Delaunay")
        if rº != sº {
          // Yes
          // Add vertex rº on the edge f
          // This will yield
          let fº = addVertex(outside: f, vertex: rº)
          showVoronoi(label: "fo_", type: "Delaunay")
  
          // Add sº (which may just be rº)
          // eº     is the edge h  -> g
          // eº + 2 is the edge sº -> h 
          let eº = addVertex(outside: e, vertex: sº)
  
          // The values of e & f can be updated here
          e = eº + 1
          f = fº + 2
  
        } else { // rº == sº
          // Reflex hubs look like this
          //
          //     g ----------- X
          //      \  eº+2     / \     
          //       \   eº+1 /e+2 \  
          //        \eº   /    e+1\ 
          //         \  /     e    \
          //          rº --------- h
          //        / \     fº   /   
          //      / f+2\fº+1   /
          //    /f      \fº+2/
          //  /   f+1    \ /     
          // j----------- Y
          //
          
          // Add the external vertices
          // fº     is the edge j -> h
          // fº + 1 is the edge h -> rº
          // No reflex or colinear
          // Add vertex rº on the edge f
          // This will yield
          var replacedEdges = addVertex(on: f, vertex: rº)
          
          // New boundary edge is fº
          let fº = replacedEdges[2] // Magic number, returns [f + 1, fº + 2, fº, f]
          showVoronoi(label: "fo_", type: "Delaunay")
          
          // Add sº (which may just be rº)
          // eº     is the edge h  -> g
          // eº + 2 is the edge sº -> h
          replacedEdges = addVertex(on: e, vertex: sº)
          let eº = replacedEdges[2] // Magic number, returns [e + 1, eº + 2, eº, e]
          
          // Need to link rº and sº since they're equal
          // Link edges - reflex corner or colinear
          link(fº, e)
          
          // This changes the hull!
          loopEdges(delete: e)
          loopEdges(delete: fº)
          
          // e needs to be updated
          e = eº
        }

        // The triangulation is no longer locally Delaunay
        var trackedEdges = Array<Int>([e, f])
        
        // Shorthand
        func updateBoundary() {
          e = trackedEdges[0]
          f = trackedEdges[1]
        }
        
        // Update the triangulation
        showVoronoi(label: foldLabel, type: "Delaunay")
        if showMe > 1 { 
          print("Flip external vertex \(rº)")
        }
        if flipVertex(onEdge: f, trackBoundaries: &trackedEdges) { updateBoundary() }
        showVoronoi(label: foldLabel, type: "Delaunay")
        
        // Convex corner
        if rº != sº {
          if showMe > 1 { 
            print("Flip external vertex \(sº)")
          }

          if flipVertex(onEdge: e, trackBoundaries: &trackedEdges) { updateBoundary() }
          showVoronoi(label: foldLabel, type: "Delaunay")
        }

        // Add the internal vertices
        // Add pº near the incoming edge
        var shell = addVertex(near: e, vertex: pº)!
        var a = shell.first!
        var a0 = 3 * (a / 3)
        
        // This can break f
        // e isn't used again but it won't be changed by an internal vertex
        // fix the boundary
        func update(edge f:Int, on shell:Array<Int>) -> Int {
          // Either added inside a triangle or on an edge
          if 3 == shell.count {
            // Added inside
            // Was f disturbed?
            for i in 1...2 {
              if f == a0 + (a + i) % 3 {
                trackedEdges[1] = shell[i]
                return trackedEdges[1]
              }
            }
          } else {
            // Added on the edge
            if f == a0 + (a + 2) % 3 { // a2
              trackedEdges[1] = shell[1] // d2
              return trackedEdges[1]
            } else {
              let b = halfEdges[a]
              if b > BoundaryEdge {
                let b0 = 3 * (b / 3)
                if f == b0 + (b + 1) % 3 { // b1
                  trackedEdges[1] = shell[2] // c1
                  return trackedEdges[1]
                }
              } // End of not a boundary
            } // End of not disturbed on shell
          }
          
          return trackedEdges[1]
        }
        
        // Update f
        f = update(edge: trackedEdges[1], on: shell)
        
        // Now we can update the mesh
        var a2 = a0 + (a + 2) % 3
        if showMe > 1 { 
          print("Flip internal vertex \(pº)")
        }
        if flipVertex(onEdge: a2, trackBoundaries: &trackedEdges) {  updateBoundary() }    
        showVoronoi(label: foldLabel, type: "Delaunay")
        
        // Add qº near the outgoing edge
        if qº != pº {
          // This is a reflex vertex
          // Add this near the outgoing spoke f
          // By construction should never duplicate a point
          shell = addVertex(near: f, vertex: qº)!
          a = shell.first!
          a0 = 3 * (a / 3)
          
          // Update f
          f = update(edge: trackedEdges[1], on: shell)
          
          // Now proceed
          a2 = a0 + (a + 2) % 3
          if showMe > 1 { 
            print("Flip internal vertex \(qº)")
          }
          if flipVertex(onEdge: a2, trackBoundaries: &trackedEdges) { updateBoundary() }        
          showVoronoi(label: foldLabel, type: "Delaunay")
        }

        // Anticlockwise search for an edge connected to the hub h
        var s = f, s2 = -1
        repeat {
          // Go anticlockwise around the hub h
          let  s0 = 3 * (s / 3)
          s2 = s0 + (s + 2) % 3
          if vertices[s2] == list[hIndex] {
            break
          }
          
          // Next edge
          s = halfEdges[s2]
        } while s > BoundaryEdge && s != f
        
        // Record an edge linked to the hub vertex we wish to remove
        // Edges are not stable in the long term
        if vertices[s2] != list[hIndex] {
          throw triangulationError.initError( "Cannot find edge connected to hub \(list[hIndex])")
        }
        
        // The hub edge
        let hubEdge = s2
    
        
        // Find first loop edge which will survive -
        // Start at rº; which is recorded on f
        a = Triangulation.hullNext[f]!
        repeat {
          let a0 = 3 * (a / 3)
          let a2 = a0 + (a + 2) % 3
          
          // Have we found a loop edge which is not connected to the hub?
          if list[hIndex] != vertices[a2] { break }
          a = Triangulation.hullNext[a]!
        } while a != f
        if a == f {
          throw triangulationError.initError("Cannot find any edge not connected to hub \(list[hIndex])")
        }
        
        // remove this vertex
        try! removeVertex(from: hubEdge)
        
        if showMe > 0 {
          countHubs += 1
          if countHubs > hubThreshold {
            hubThreshold += hubStep
            print("Processed \(countHubs) hubs : \(100.0 * Double(countHubs) / numberHubs)%")
          }          
        }
        
        // Get next edge
        e = a
        repeat {
          f = e
          e = Triangulation.hullPrev[f]!
          if e == a {
            throw triangulationError.initError("Cannot find edge connected to rº \(rº)")
          }
        } while vertices[e] != rº
        
        // Next list
        if vertices[f] != loopStop {
          list = loopVertices[iº]!
        }
      } while vertices[f] != loopStop
    } // Finished all loops
    
    // Fifth Pass
    // Balance the remaining loop edges by adding image vertices
    if showMe > 0 {
      print("Pass 5")
      print("\tAdd Image Vertices")
    }
    try! addImages()
  }
  
  // Flip around a vertex
  mutating func flipVertex(onEdge e:Int, trackBoundaries boundaryList:inout Array<Int>) -> Bool {
    // Flip vertex
    var flipped = false
    // Maintains a stack of edges for flipping
    // Differs from flipEdges in that instead
    // of considering a triangle at a time it considers the
    // vertex attached to the flipped edge and adds that to the stack
    
    //                     *       *
    //                      \     /
    //                       \   /
    //                        \ /  e
    //                * ------ v ----- *
    //                        / \
    //                       /   \
    //                      /     \
    //                     *       *
    // Establish the stopEdge
    var stopEdge = EmptyEdge
    
    // Now get the vertex
    var v = vertices[e]
    
    // The initial edge
    var a = e
    
    // Now set the vertex stack
    var vertexStack = Dictionary<Int, Int>()

    // eachVertex : repeat {
    eachVertex: while EmptyVertex != v {
      var a2:Int
      
      eachEdge: repeat {
        let a0 = 3 * (a/3)
        let a1 = a0 + (a + 1) % 3
        a2 = a0 + (a + 2) % 3
        
        // Start by trying to flip edge a
        if let flipList = flipLeft(edge: a, unconditional: false) {
          /* flipped a
           *
           *            i                   i
           *           /  \               /| |\
           *        a2/    \a1        a2 / | | \b2
           *         /      \           /  | |  \
           *        /___a____\         /   | |   \
           *       v__________w  =>   v  a1| |b   w
           *        \  b2    /         \   | |   /
           *         \      /           \  | |  /
           *        b \    /b1         a \ | | /b1
           *           \  /               \| |/
           *            p                   p
           *
           *
           */
          flipped = true
          
          // Update the vertex stacks
          // Only one entry is needed for each vertex
          let p  = vertices[a1], i = vertices[a2]
          let b  = halfEdges[a1]
          let b0 = 3 * (b / 3)
          let b2 = b0 + (b + 2) % 3
          let w = vertices[b2]
          if showMe > 1 { print("\tFlipped edge connecting \(v) => \(w)") }
          
          // We are processing (v) so can be ignored
          // If w is in the stack replace the edge
          if nil != vertexStack[w] {
            vertexStack[w] = b2
          }
          // Now add/replace entries for (i) and (p)
          vertexStack[i] = b
          vertexStack[p] = a1
          
          // Keep the boundary edge eº valid
          // a2 and b1 are not disturbed,
          // a and b2 cannot be boundary edges
          if !flipList.isEmpty {
            for (index, tracked) in boundaryList.enumerated() {
              // Was this edge replaced
              if let newEdge = flipList[tracked] {
                boundaryList[index] = newEdge
              }
            }
          }
          
          // If stopEdge was changed update it
          if b == stopEdge {
            stopEdge = a
          }
          // Repeat with edge a
          continue eachEdge
        } else if let flipList = flipLeft(edge: a1, unconditional: false) {
          // If the spoke did not flip try the shell
          // This creates an extra spoke if it works
          /* flip an edge a
           *
           *            i                     i
           *          /| |\                  /  \
           *      a2 / | | \b1            b2/    \b1
           *        /  | |  \              /      \
           *       /   | |   \    flip    /___b____\
           *      v  a1| |b2   w     =>  v__________w
           *       \   | |   /            \   a2   /
           *        \  | |  /              \      /
           *       a \ | | /b              a\    /a1
           *          \| |/                  \  /
           *            p                     p
           *
           */
          flipped = true

          // Update the vertex stacks
          // Only one entry is needed for each vertex
          let w = vertices[a2]
          if showMe > 1 { print("\tFlipped edge connecting \(v) => \(w)") }
          
          // We are processing (v) so can be ignored
          // Now add/replace entry for (w)
          vertexStack[w] = a2
          
          // Keep the boundary edge eº valid
          // a and  b1 are not disturbed,
          // a1 and b2 cannot be boundary edges
          if !flipList.isEmpty {
            for (index, tracked) in boundaryList.enumerated() {
              // Was this edge replaced
              if let newEdge = flipList[tracked] {
                boundaryList[index] = newEdge
              }
            }
          }
        } else  if stopEdge == EmptyEdge {
          // Edge a is made the stopEdge
          stopEdge = a
        }
        
        // Next edge if no flip
        a = halfEdges[a2]
        
        // Test for stop edge here - to prevent premature break
        if a == stopEdge {
          break
        }
        
        // Stop at the stopEdge
      } while a > BoundaryEdge
      
      if showMe > 1  && flipped { print("Flip Vertex \(v)") }
      
      // Now get another vertex
      // Query : Does order matter?
      if let first = vertexStack.first {
        // Establish the stopEdge
        stopEdge = EmptyEdge
        
        // And the edge
        v = first.key
        a = vertexStack.removeValue(forKey: v)!
      } else {
        v = EmptyVertex
      }
    }
    
    // Any processing?
    return flipped
  }
  
  // This adds the new vertices along each active edge (from one rº => sº)
  mutating func addImages() throws {
    showMe -= 1
    // Get each hull loop
    for (aº, bº) in Triangulation.loopHooks {
      // Find the edge joining the vertices aº & bº
      let loopStart = findLoopEdge(from: aº, to: bº)
            
      // So the stop vertex is aº
      let loopStop = aº
      
      // Get a boundaryedge from this loop
      //
      //      qº              pª (j)
      //      |               |
      //   h •|•••••••••••••••|•j
      //      |       eº      |
      //      rº--------------sª (j)
      //
      
      // Since loop hook always connected to convex hub this will work
      var vertexList =  Triangulation.loopVertices[bº]
        
      // Unpack the vertices
      var qª = vertexList![qIndex], pª = vertexList![pIndex],
          rª = vertexList![rIndex], sª = vertexList![sIndex]
      var j  = vertexList![hIndex]
      
      // Set the incoming edge to be from sª to rª
      var eº = loopStart
      
      // The logic here is based on edge segments rather than vertices
      roundLoop: repeat {
        // From the incoming edge generate the outgoing one
        // Can't here assume that the hub is convex;
        // can be reflex or colinear
        if rª != sª {
          // Convex
          // Move forward
          eº = Triangulation.hullNext[eº]!
        }
        
        // We can set h
        let h = j
        
        // And the left hand sentinels
        let rº = rª, qº = qª

        // Find the right hand sentinels
        var fº = Triangulation.hullNext[eº]!

        // If this is a convex hub this will be empty
        vertexList = Triangulation.loopVertices[vertices[fº]]
        while nil == vertexList {
          fº = Triangulation.hullNext[fº]!
          vertexList = Triangulation.loopVertices[vertices[fº]]
        }
        
        // The vertices associated with the previous entry
        j  = vertexList![hIndex]
        
        // Get the sentinels
        // Reflex and Convex corners different
        qª = vertexList![qIndex]; pª = vertexList![pIndex]
        rª = vertexList![rIndex]; sª = vertexList![sIndex]
        
        // These are fixed for the time being
        let hx = Triangulation.coords[2 * h], hy = Triangulation.coords[2 * h + 1]
        let jx = Triangulation.coords[2 * j], jy = Triangulation.coords[2 * j + 1]
        
        // The vertices that define  the zip line
        var q = qº, p = pª
        var r = rº, s = sª
        
        // Not all these special cases should occur
        showVoronoi(label: "be_")
        var mirrorVertices: Dictionary<Int, Int> = [s : p]
        repeat {
          // Balance one vertex at a time
          
          //    Five Types (ignoring the co-linear special cases!)
          //
          // No diagonal: zero flip  Right Diagonal: Single Flip  Left Diagonal: Single Flip
          //             v                      v----- pª                     v----- pª
          //            / \                     |    / |                      | \    |
          //    A      /   \              B     |   /  |              C       |  \   |
          //          /     \                   |  /   |                      |   \  |
          //         /       \                  | /    |                      |    \ |
          //       rº____vº___sª                rº_____sª                     rº____sª
          //
          //             Double Flip
          //        Left Diagonal  Right Diagonal
          //
          //            v              v
          //           / \            / \
          //          /   \          /   \
          //         /     \        /     \
          //        qº______pª     qº______pª
          //        | \     |      |     / |
          //        |  \    |      |    /  |
          // E & D  |   \   |      |   /   |
          //        |    \  |      |  /    |
          //        |     \ |      | /     |
          //        rº______sª     rº______sª
          //
          //
          
          var v = EmptyVertex
          var unbalanced = false
          var a = eº
                              
          // Search Anti-clockwise around rº
          // Can find
          // (A) Zero Flip
          // or (B) Single Flip (Right)
          //
          
          // Set diagonal properties
          var d = 0
          var rightDiagonal = EmptyEdge
          var leftDiagonal = EmptyEdge

          repeat {
            let a0 = 3 * (a/3)
            let a2 = a0 + (a + 2) % 3
            
            // Candidate vertex
            v = vertices[a2]
            let ha2 = halfEdges[a2]
            if ha2 <= BoundaryEdge {
              throw triangulationError.initError("Found an unexpected boundary near edge \(eº) anticlockwise around \(r)")
            }
            // Done when we reach q or (when q == r) a boundary
            if v == q {
              rightDiagonal = ha2
              v = EmptyVertex
              break
            }
            
            // Encroaching if not p
            if v != p {
              // Is this vertex unbalanced?
              unbalanced = insideCircumCircle(vertex: v, p, q, r)
              
              // Handle type A
              if rightDiagonal == EmptyEdge {
                if !unbalanced {
                  
                  // Type A vertex - Balanced Not Encroaching
                  // Convert it to type C
                  // A->C does not change boundary
                  // The flip must be forced
                  _ = flipRight(edge: a2, unconditional: true)
                  
                  // Record the new left diagonal
                  d = 2
                  leftDiagonal = a0 + (a + 1) % 3
                  a = leftDiagonal
                  
                  // Type C Vertex
                  showVoronoi(label: "be_")
                } // else Unbalanced Type A
              } // else Type B (Balanced or Unbalanced)
              
              // At this stage
              // A (Unbalanced) or B
              // Not rightDiagonal
              break
            }
            
            // The diagonal r->p exists
            // So no encroaching vertices can be connected to sª
            rightDiagonal = ha2
            d = 1
            
            // Next edge
            a = rightDiagonal
          } while true
          
          // If no vertex found search around the right hand pivot (s)
          // Only needed if (p) was not seen
          if EmptyVertex == v && EmptyEdge == rightDiagonal {
            // Clockwise around sª
            a = eº
            repeat {
              let a0 = 3 * (a/3)
              let a1 = a0 + (a + 1) % 3
              let a2 = a0 + (a + 2) % 3
              
              // The vertex is on spoke a2; we advance through spoke a1
              v = vertices[a2]
              let ha1 = halfEdges[a1]
              if ha1 <= BoundaryEdge {
                throw triangulationError.initError("Found an unexpected boundary near edge \(eº) clockwise around \(s)")
              }
              
              // Done when we reach p
              if v == p {
                leftDiagonal = ha1
                v = EmptyVertex
                break
              }
              
              // Encroaching if not q
              if v != q {
                // Is this vertex unbalanced?
                unbalanced = insideCircumCircle(vertex: v, p, q, r)

                // At this stage
                // C
                break
              }
              
              // The diagonal q->s exists
              leftDiagonal = ha1
              d = 2
              
              // Next edge - will never be a boundary since
              // the search stops no later than p
              a = leftDiagonal
            } while true
          }
          
          // If no vertex found check the triangle connected to (q->p)
          if leftDiagonal != rightDiagonal {
            if !unbalanced {
              // A balanced type B or C vertex was found
              // All type A are unbalanced at this stage
              //
              // Convert these to type D or E
              // Reverse the offset
              if EmptyEdge != rightDiagonal {
                a = rightDiagonal
                d = 2
              } else {
                a = leftDiagonal
                d = 1
              }
              
              let a0 = 3 * (a/3)
              a = a0 + (a + d) % 3
              
              // Why is this always a flip right?
              _ = flipRight(edge: a, unconditional: false)
              
              // No candidate
              v = EmptyVertex
            } else if EmptyVertex == v {
              // Type B or C
              // Check for types D & E
              // One of the diagonals exists - since when A it is unbalanced, so v != EmptyVertex
              // Check for balanced a
              a = 1 == d ? rightDiagonal : leftDiagonal
              
              // Must be type D or E
              // Get the edge q->p
              var a0 = 3 * (a/3)
              a = a0 + (a + d) % 3
              let pq = halfEdges[a]
                
              // Get the opposite vertex
              a0 = 3 * (pq/3)
              let a2 = a0 + (pq + 2) % 3
              
              // The candidate vertex
              v = vertices[a2]
              
              // Either D or E
              
              // If it is unbalanced convert this to type B or C
              // This is done by flipping the edge pq
              // BUT don't force a flip - since this
              // can make a non-delaunay mesh
              var didFlip = false
              
              if insideCircumCircle(vertex: v, p, q, r) {
                // Left diagonal would be undisturbed
                // Right diagonal is changed by single flip
                if rightDiagonal != EmptyEdge {
                  didFlip = (nil != flipLeft(edge: pq, unconditional: false))

                  a = rightDiagonal
                } else {
                  // Use a right flip
                  didFlip = flipRight(edge: pq, unconditional: false)

                  a = leftDiagonal
                }
              }
              
              // Did it flip?
              if !didFlip {
                // Not a candidate
                // Left as D or E
                v = EmptyVertex
              }
            }
          }
          
          // At this point we have:
          //
          // A - Unbalanced Candidate
          // B - Candidate
          // C - Candidate
          // D - No Candidate
          // E - No Candidate
          
          // If a candidate exists balance it
          var didFlip = false
          if EmptyVertex != v {
            // Found a candidate; so A, B or C

            // We have either type A, B or C
            // Convert type B & C to type A
            if d > 0 {
              // We have the diagonal already
              if EmptyEdge != rightDiagonal {
                // Convert B->A, use right flip
                didFlip = flipRight(edge: a, unconditional: false)

              } else {
                // Left diagonal
                // Convert C->A with left flip does not disturb boundary
                didFlip = (nil != flipLeft(edge: a, unconditional:false))
              }
              
              // If flip failed then retreat (uncondiotonall
              if !didFlip {
                // Convert to B,C -> D, E
                if EmptyEdge != rightDiagonal {
                  a = rightDiagonal
                  d = 2
                } else {
                  a = leftDiagonal
                  d = 1
                }
                
                let a0 = 3 * (a/3)
                a = a0 + (a + d) % 3
                
                // A right flip
                _ = flipRight(edge: a, unconditional: true)
                
                // No candidate
                v = EmptyVertex
              }
              
              showVoronoi(label: "be_")
            } else {
              // Class A - treat as "did flip"
              didFlip = true
            }
          }
          
          if (didFlip) {
            // Add vertex v *on* eº
            //
            //        v
            //       /|\
            //      / | \
            //     /  |  \
            //    r---vº--s
            //        eº
            //
            
            // Reflect this vertex in the line h->j
            let vx = Triangulation.coords[2 * v], vy = Triangulation.coords[2 * v + 1]
            
            // Reflect this vertex
            let vº = reflect(vx, vy, hx, hy, jx, jy)
            
            // The local edge code
            Triangulation.properties.append(Triangulation.properties[v])
            growTriangulation(to: Triangulation.pointCount)
            
            // fº + 2 is the edge v  -> rº
            // fº     is the edge rº -> vº
            // fº + 1 is the edge vº -> x
            let replacedEdges = addVertex(on: eº, vertex: vº)
            showVoronoi(label: "be_")

            // New boundary edge is fº
            let fº = replacedEdges[2] // Magic number, returns [eº + 1, fº + 2, fº, eº]

            // Update the edge
            eº = fº
            
            // Need to flip edges connected to v
            // Because various flips were made to get it to
            // the sought configuration
            // v is found on replacedEdges[eº+2] = fº+2
            var trackedEdges = Array<Int>([eº])
            if flipVertex(onEdge: fº + 2, trackBoundaries: &trackedEdges) {
              eº = trackedEdges.first!
            }
            
            // A new vertex - q & r are undisturbed
            // But save the old p
            mirrorVertices[s] = p
            
            p = v
            s = vº
          } else {
            // No candidate; D or  E
            // Move one edge segment along
            q = p
            
            // Next edge
            eº = Triangulation.hullNext[eº]!
            
            //
            r = s
            s = vertices[Triangulation.hullNext[eº]!]
            p = mirrorVertices[s, default:pª]

          }
        } while q != p
        
        // Is this the stop vertex
      } while loopStop != sª
    } // End of all loops
  }
  
  // Remove a vertex from the triangulation
  mutating func removeVertex(from startEdge:Int) throws {
    // Nested functions
    func unlink(_ a:Int) -> Int {
      // Unlinking always sets the default boundary code
      // if it was not already a boundary code
      let ha = halfEdges[a]
      if ha > BoundaryEdge {
        // A is labelled as empty
        halfEdges[ha] = EmptyEdge
      }
            
      // Return opposing half edge
      return ha
    }
    
    // Remove any vertex from the triangulation
    let hVertex = vertices[startEdge]

    // The initial edge
    var a = startEdge

    // Establish the stopEdge (might be a boundary code)
    // We have to note it because triangles will be removed
    let stopEdge = halfEdges[startEdge]
        
    // Build a list of vertices
    var shellVertices = Array<Int>()
    
    // And a list of the corresponding edges
    var shellEdges = Array<Int>()
    
    // Go anticlockwise starting at p
    //
    //                       t -----s
    //                      /| d1  /|
    //                     / |    / |
    //                    / e|d2 /  |
    //                   /e1 | d/c2 |
    //                  /    | /  c1|
    //                 /  e2 |/  c  |
    //                u------h------r
    //                      / \ b2  |
    //                     /   \  b1|
    //                    /   a2\   |
    //                   /a      \b |
    //                  /         \ |
    //                 /    a1     \|
    //                p ----------- q
    //
    //
    //    Need to update the triangulation loops
    //
    // Two possibilities; h is on the loop (a convex corner)
    //                    h is not on the loop (reflex or colinear)
    if stopEdge <= BoundaryEdge {
      shellEdges.append(startEdge)
    }
    
    var a2:Int
    vertexLoop: repeat {
      let a0 = 3 * (a / 3)
      let a1 = a0 + (a + 1) % 3
      a2 = a0 + (a + 2) % 3
      
      // Pivot around h 
      let ha2 = halfEdges[a2]
      let v = vertices[a2]

      // Remove the triangle
      // Unlink spokes
      _ = unlink(a)
      _ = unlink(a2)

      // Save removed triangle id
      emptyTriangles.insert(a0)
      
      // Save shell edge and vertex
      shellVertices.append(vertices[a1])
      shellEdges.append(a1)

      // Is this the stopEdge - can't rely on the halfEdge
      // still being connected
      if stopEdge == ha2 || ha2 <= BoundaryEdge {
        // If this is a boundary edge the
        // last vertex u has not been added
        // Save it
        
        if stopEdge <= BoundaryEdge {
          shellVertices.append(v)
          shellEdges.append(a2)
        }
        
        // All done
        break vertexLoop
      }
      
      // Get the next triangle
      a = ha2
      
      // Stop here if edge wasn't disconnected (it will have been)
    } while a2 != stopEdge


    if showMe > 1 {
      print("Remove Vertex => \(vertices[startEdge])")
      print("\tShell Edges    => \(shellEdges)")
      print("\tShell Vertices => \(shellVertices)")
    }
    
    // Logging
    showVoronoi(label: "rv_\(hVertex)_", type:"Delaunay")
    
    // Get each outer edge
    // They will either be *boundary* edges or *internal* edges
    // We need to add all the internal edges to the loops, connecting them
    // up in the right way; and also remove any boundary edges

    // How the edges are to be processed - essentially edges are turned inside out
    //
    //     . d          f .            .\ d         f /.
    //   g  .            . i           g.\           /.i
    //       . ...e.... .                .\         /.
    //       .*––––––––*.           ==>   .*       *.
    //       .|    h   |.                 .|       |.       n was Triangulation.hullNext[e]
    //      p.|        |.n               p.|       |.n      p was Triangulation.hullPrev[e]
    //
    //
    //        ....
    //        ____ boundary (inside is up)
    //
    //        .... internal
    //
    //   First pass skips boundaries
    //
    //   Possible pairs are (B = Boundary, I = Internal, F = Flipped Boundary (Former Internal)
    //   1) B*B
    //   2) B*I
    //   3) F*I
    //   4) F*B
    //   5) I*I (First edge only)
    //   6) I*B (First edge only)
    //
    //  Type 6 can only occur if stopEdge was a BoundaryEdge
    //  Rotating the shell will put B as the last edge
    //  So possible initial pair will be
    //   i) B*B (== 1)
    //  ii) B*I (== 2)
    var bLast = false

    // Initialize
    // Get the previous edge (will be type B or I)
    var d = shellEdges.last!
    if halfEdges[shellEdges.first!] <= BoundaryEdge {
      // Put this B at the end (so I*B is impossible)
      d = shellEdges.removeFirst()
      shellEdges.append(d)
    }
    
    // The corresponding edge
    var g = halfEdges[d]
                
    // Now process each shell edge in turn (anti-clockwise order)
    for (j, e) in shellEdges.enumerated() {
      // Previous edge pair is (d, g)
      // This edge pair is (e, h)
      
      // We have not unlinked these edges yet
      let h = unlink(e)
      
      // Possible options
      // I*I, B*I, F*I
      if h > BoundaryEdge {
        // Type F*I 
        if g <= BoundaryEdge {
          // Type B*I
          let n = Triangulation.hullNext[d]!
          Triangulation.hullNext[h] = n
          Triangulation.hullPrev[n] = h
          
          // Can delete d if not first edge
          if j > 0 {
            Triangulation.hullNext[d] = nil
            Triangulation.hullPrev[d] = nil
          } else {
            bLast = true
          }
        } else if d <= BoundaryEdge || 0 == j {
          // Type F*I or I*I
          Triangulation.hullNext[h] = g
          Triangulation.hullPrev[g] = h          
        } else  {
          // I*I can only occur on first call
          throw triangulationError.initError("Found unexpected edge pair status")
        }
        
        // Set the edge code
        halfEdges[h] = BoundaryEdge

        // Update edge
        d = BoundaryEdge
      } else {
        // Type *B
        if d <= BoundaryEdge {
          // Type F*B
          let p = Triangulation.hullPrev[e]!
          Triangulation.hullPrev[g] = p
          Triangulation.hullNext[p] = g
        } else if g <= BoundaryEdge {
          // Type B*B
          if 0 != j {
            // Nothing to do for B*B
            // Can delete d (unless first call)
            Triangulation.hullNext[d] = nil
            Triangulation.hullPrev[d] = nil
          } else {
            bLast = true
          }
        } else {
          // This is I*B
          // Not B*B or F*B
          throw triangulationError.initError("Found unexpected edge pair status")
        }

        // Update edge
        d = e

        // Logging
        showVoronoi(label: "rv_\(hVertex)_", type:"Delaunay")
      }
      
      // Update edges
      g = h
    }
    
    // Delete boundary
    if bLast {
      d = shellEdges.last!
      Triangulation.hullNext[d] = nil
      Triangulation.hullPrev[d] = nil
    }
    
    // Logging
    showVoronoi(label: "rv_\(hVertex)_", type:"Delaunay")
    
    // Now use the shell of vertices to define a zone - use full machinery to reattach
    // At the moment this doesn't add any extra points but does
    // add vertices and triangles
    // treat it as if it did add points
    let localZone = Zone(given: shellVertices)

    // Consider reusing empty triangles to avoid this step
    Triangulation.pointCount += localZone.polygonCount
    growTriangulation(to: Triangulation.pointCount)

    // Next triangulate this zone
    for z in localZone.convexZones {
      // Track this
      
      // Add the zone to the triangulation
      let convexHullNext = try! addVertices(indices: z.zoneIndices)
      
      // Join convexZones together
      //
      try! join(loop: convexHullNext)
      showVoronoi(label: "rv_\(hVertex)_", type: "Delaunay")
    } // End of each convex zone
  }
  
  // Flip one edge
  // Right flip
  mutating func flipRight(edge a:Int, unconditional forceFlip:Bool = false) -> Bool {
    /* if the pair of triangles doesn't satisfy the Delaunay condition flip it
     *
     *    Input edge is a
     *
     *
     *            h                     h
     *          /| |\                  /  \
     *      a1 / | | \b2            a1/    \a
     *        /  | |  \              /      \
     *       /   | |   \    flip    /___a2___\
     *      i   a| |b  p     =>    i__________p
     *       \   | |   /            \   b2   /
     *        \  | |  /              \      /
     *       a2\ | | /b1             b\    /b1
     *          \| |/                  \  /
     *            j                     j
     *
     *
     */
    let b:Int = halfEdges[a]
    if b <= BoundaryEdge { return false }

    // Get the edges associated with each triangle
    // Triangle a - need a2 for special case
    let a0 = a - a % 3
    let a2 = a0 + (a + 2) % 3

    // Finish Triangle a & Triangle b
    let b0 = b - b % 3
    let b2 = b0 + (b + 2) % 3
    
    // Now the swapped points
    let i = vertices[a2]
    let p = vertices[b2]
  
    if !forceFlip {
      // Need to check if this flip maintains the locally Delaunay property
      //
      let h = vertices[b]
      let j = vertices[a]

      // We need to flip the edges a-b if the point p is inside the circum-circle
      // of the triangle [h, i, j]
      
      // Use Shewchuk's robust version of inCircumCircle
      let illegal = inCircumCircle(
        ax: Triangulation.coords[2 * h], ay: Triangulation.coords[2 * h + 1],
        bx: Triangulation.coords[2 * i], by: Triangulation.coords[2 * i + 1],
        cx: Triangulation.coords[2 * j], cy: Triangulation.coords[2 * j + 1],
        px: Triangulation.coords[2 * p], py: Triangulation.coords[2 * p + 1])
      
      // If not illegal all done
      if !illegal { return false }
    }
    
    // Proceed with the flip
    vertices[a] = p
    vertices[b] = i
        
    // Link the half edges
    //   halfEdges[a]  <=> halfEdges[b2]
    //   halfEdges[a2] <=> halfEdges[b2]
    //   halfEdges[b]  <=> halfEdges[a2]
    let ha2 = halfEdges[a2], hb2 = halfEdges[b2]
    link(a, hb2)
    link(b, halfEdges[a2])
    link(a2, b2)
        
    // This can perturb boundary at b2 => a and a2 => b
    if hb2 <= BoundaryEdge {
      loopEdges(edge: a, replaces: b2)
    }
    if ha2 <= BoundaryEdge {
      loopEdges(edge: b, replaces: a2)
    }
    
    // Return that it flipped
    return true
  }
  
  // Left flip
  mutating func flipLeft(edge a:Int, unconditional forceFlip:Bool = false) -> Dictionary<Int, Int>? {
    /* flip an edge a
     *
     *    Input edge is a - need not to disturb neighbouring spoke b2
     *
     *
     *            h                     h
     *          /| |\                  /  \
     *      a1 / | | \b2             b/    \b2
     *        /  | |  \              /      \
     *       /   | |   \    flip    /___b1___\
     *      i   a| |b  p     =>    i__________p
     *       \   | |   /            \   a1   /
     *        \  | |  /              \      /
     *       a2\ | | /b1            a2\    /a
     *          \| |/                  \  /
     *            j                     j
     *
     *
     */
    var flipList = Dictionary<Int, Int>()
    
    let b = halfEdges[a]
    if b <= BoundaryEdge { return nil }
    
    // Get the edges associated with each triangle
    // Triangle a
    let a0 = a - a % 3
    let a1 = a0 + (a + 1) % 3
    let a2 = a0 + (a + 2) % 3
    
    // Finish Triangle a & Triangle b
    let b0 = b - b % 3
    let b1 = b0 + (b + 1) % 3
    let b2 = b0 + (b + 2) % 3
    
    // Now the swapped points
    let i = vertices[a2]
    let p = vertices[b2]
    
    if !forceFlip {
      // Need to check if this flip maintains the locally Delaunay property
      //
      let h = vertices[b]
      let j = vertices[a]
      
      // We need to flip the edges a-b if the point p is inside the circum-circle
      // of the triangle [h, i, j]
      
      // Use Shewchuk's robust version of inCircumCircle
      let illegal = inCircumCircle(
        ax: Triangulation.coords[2 * h], ay: Triangulation.coords[2 * h + 1],
        bx: Triangulation.coords[2 * i], by: Triangulation.coords[2 * i + 1],
        cx: Triangulation.coords[2 * j], cy: Triangulation.coords[2 * j + 1],
        px: Triangulation.coords[2 * p], py: Triangulation.coords[2 * p + 1])

      // If not illegal all done
      if !illegal { return nil }
    }
    
    // Proceed with the flip
    vertices[a1] = p
    vertices[b1] = i
    
    // Link the half edges
    //   halfEdges[a]  <=> halfEdges[b1]
    //   halfEdges[b]  <=> halfEdges[a1]
    //   halfEdges[a1] <=> halfEdges[b1]
    let ha1 = halfEdges[a1], hb1 = halfEdges[b1]
    link(a, hb1)
    link(b, ha1)
    link(a1, b1)
    
    // This can perturb boundary at b1 => a and a1 => b
    if hb1 <= BoundaryEdge {
      loopEdges(edge: a, replaces: b1)
      
      // Means b1 replaced by a
      flipList[b1] = a
    }
    if ha1 <= BoundaryEdge {
      loopEdges(edge: b, replaces: a1)
      
      // Means a1 replaced by b
      flipList[a1] = b
    }
    
    // Return that it flipped
    return flipList
  }

  // is a vertex in a loop
  func loopEdges(around vertex:Int) -> Array<Int> {
    return Triangulation.hullNext.compactMap { vertex == vertices[$0.key] ? $0.key : nil }
  }
    
  // find an edge joining two adjacent and ordered loop vertices
  func findLoopEdge(from s:Int, to r:Int) -> Int {
    if s != r {
      return Triangulation.hullNext.first(where: {s == vertices[$0.key] && r == vertices[$0.value]})!.key
    }
    
    // Just get the first loop edge connected to s
    return Triangulation.hullNext.first(where: {s == vertices[$0.key]})!.key
  }
  
  // Delete edge e from the boundary
  func loopEdges(delete e:Int) {
    // The edge e is deleted
    //       d   |   e   |   f
    //     ----- X       X ----->
    let d = Triangulation.hullPrev[e]!
    let f = Triangulation.hullNext[e]!
    Triangulation.hullPrev[f] = d;     Triangulation.hullNext[d] = f
    Triangulation.hullPrev[e] = nil;   Triangulation.hullNext[e] = nil
  }
  
  // Insert edge a into the boundary after the edge e
  func loopEdges(insert a:Int, after e:Int) {
    // The edge a is inserted after e
    //       e   |   a    |   f
    //     ----- V ------ V ----->
    
    let f = Triangulation.hullNext[e]!
    Triangulation.hullPrev[a] = e;     Triangulation.hullNext[a] = f
    Triangulation.hullPrev[f] = a;     Triangulation.hullNext[e] = a
  }
  
  // Insert edge a into the boundary before the edge e
  func loopEdges(insert a:Int, before e:Int) {
    // The edge a is inserted before e
    //       d   |   a    |   e
    //     ----- V ------ V ----->
    
    let d = Triangulation.hullPrev[e]!
    Triangulation.hullPrev[a] = d;     Triangulation.hullNext[a] = e
    Triangulation.hullPrev[e] = a;     Triangulation.hullNext[d] = a
  }
  
  // Insert edge a into the boundary replacing the edge e
  func loopEdges(edge a:Int, replaces e:Int) {
    // The edge a replaces e
    //       d   |   a    |   f
    //     ----- V ------ V ----->
    let d = Triangulation.hullPrev[e]!
    let f = Triangulation.hullNext[e]!
    Triangulation.hullPrev[a] = d;     Triangulation.hullNext[a] = f
    Triangulation.hullPrev[f] = a;     Triangulation.hullNext[d] = a
    Triangulation.hullPrev[e] = nil;   Triangulation.hullNext[e] = nil
  }
  
  // Simply add a vertex ouside the boundary edge a
  mutating func addVertex(outside a:Int, vertex i:Int) -> Int {
    //
    // add the triangle from the point          i
    // This is actually just [e, i, q]         / \
    //                                        q---e
    //                                          a
    // The linked edges are BoundaryEdge (outside the hull) for (e->i) and (i->q) and hullTri[e] for (q->e)
    // The edge a is on the existing triangle a (e->q)
    // In the new triangle  the edges are (q->e) (e->i) (i->q)
  
    // Add the new vertex i to the triangulation
    // Other initial edges
    // Compute next edge
    let a0 = a - a % 3
    let a1 = a0 + (a + 1) % 3
    
    // Save the vertices
    let e = vertices[a]
    let q = vertices[a1]

    // New triangle
    let t = addTriangle(q, e, i, a, BoundaryEdge, BoundaryEdge)
    
    // a (e->q) was a boundary edge; and is replaced by t+1 (e->i) and t+2 (i->q)
    //
    //                i
    //               / \
    //        t+2   /   \ t+1
    //             /  t  \
    //      ----- q <----- e -------
    //       s        a        u
    //
    loopEdges(insert: t + 1, before: a)
    loopEdges(edge: t + 2, replaces: a)
    
    return t
  }
  
  mutating func addVertex(inside a:Int, vertex p:Int) -> Array<Int> {
    // Add the new vertex p to the triangulation
    // Other initial edges
    // Compute next edge
    let a0 = a - a % 3
    let a1 = a0 + (a + 1) % 3
    let a2 = a0 + (a + 2) % 3
    
    // Save old half edge links
    let ha1 = halfEdges[a1], ha2 = halfEdges[a2]
    
    // Save the vertices
    let q = vertices[a]
    let i = vertices[a2]
    let n = vertices[a1]
    
    // Insert vertex p on edge a
    vertices[a2] = p
    
    //
    //
    //                            i
    //                           /|\
    //                          / | \
    //                         /  |  \
    //                        /   |   \
    //                       /    |    \
    //                      /  c2 | b1  \
    //                     /      |      \
    //                    /       |       \
    //              ha2  / c      |      b \  ha1
    //                  /        / \        \
    //                 /       /  p  \       \
    //                /   c1 /         \ b2   \
    //               /     /             \     \
    //              /    / a2           a1 \    \
    //             /   /                     \   \
    //            /  /                         \  \
    //           / /              a              \ \
    //         q ----------------------------------  n
    //                           ha
    //
    //     Normal case - create new triangles b and c and reconnect a as shown
    //
    
    // Now create new triangles B and C
    // Don't assume c = b + 3 in case deleted triangles are reused
    let b = addTriangle(n, i, p, ha1, UnusedEdge, a1)
    let c = addTriangle(i, q, p, ha2, a2, b + 1)
    
    // The replaced outer edges are just b & c
    if ha1 <= BoundaryEdge {
      loopEdges(edge: b, replaces: a1)
    }
    if ha2 <= BoundaryEdge {
      loopEdges(edge: c, replaces: a2)
    }
    
    return [a, b, c]
  }
  
  // Add a vertex exactly on an edge - only happens with degenerate input in practice
  mutating func addVertex(on a:Int, vertex p:Int) -> Array<Int> {

    // Collinear - two cases; input vertex is internal or input vertex is a new hull vertex
    //   New point added *exactly* on edge joining vertices q->n; separating triangles A and B
    //
    //
    //                            i
    //                           /.\
    //                          / . \
    //                         /  .  \
    //                        /   .   \
    //                 ha2   /    .    \ ha1
    //                      /     .     \
    //                     / d2   .      \
    //                    /    d1 . a2    \
    //                   /        .     a1 \
    //                  /   d     .  a      \
    //                 q -------- p -------- n
    //                  \   c     .  b      /
    //                   \        .     b2 /
    //                    \ c1    .       /
    //                     \   c2 . b1   /
    //                      \     .     /  hb2
    //                hb1    \    .    /
    //                        \   .   /
    //                         \  .  /
    //                          \ . /
    //                            r
    //
    //     Second case is the same except edge a [q->n] is a boundary edge
    //     In this case there is no triangle B
    // Add the new vertex p to the triangulation
    
    // Initial edges
    let a0 = a - a % 3
    let a1 = a0 + (a + 1) % 3
    let a2 = a0 + (a + 2) % 3
    
    // Save old half edge links
    let ha2 = halfEdges[a2]
    
    // Save the vertices
    let q = vertices[a]
    let i = vertices[a2]
    
    // Insert vertex p on edge a
    vertices[a] = p
    
    //  Create new triangle D
    let d = addTriangle(q, p, i, BoundaryEdge, a2, ha2)
    
    // Outer edges
    if ha2 <= BoundaryEdge {
      loopEdges(edge: d + 2, replaces: a2)
    }

    // Does triangle b exist?
    let b = halfEdges[a]
    if b > BoundaryEdge {
      // Yes
      let b0 = b - b % 3
      let b1 = b0 + (b + 1) % 3
      let b2 = b0 + (b + 2) % 3
  
      // Save old half edge links
      let hb1 = halfEdges[b1]
      
      // Save the new vertex
      let r = vertices[b2]
      
      // Insert vertex p on edge b1
      vertices[b1] = p
      
      //  Create new triangle c
      let c = addTriangle(p, q, r, d, hb1, b1)
      
      // Extra Outer edges
      if hb1 <= BoundaryEdge {
        loopEdges(edge: c + 1, replaces: b1)
      }
      
      // Return outer edges
      return [a1, d + 2, c + 1, b2]
    } else {
      // a (q->n) was a boundary edge; and is replaced by
      // d (q->p) and a (p->n)
      loopEdges(insert: d, before: a)
      
      // Return outer edges
      return [a1, d + 2, d, a]
    }
    
  }
  
  // General add a vertex given only an edge hint
  mutating func addVertex(near hint:Int, vertex p:Int) -> Array<Int>? {
    var onEdge = EmptyEdge
    
    // Find a triangle enclosing a vertex p
    func findTriangle(edge hint:Int, vertex p:Int) throws -> Int {
      // Need to find which triangle the vertex lies within
      // Keep a list of all searched triangles
      var searchedTriangles = Set<Int>()
      var searchStack = Array<Int>()
      var a = hint
      
      // On first entry orientation of a with respect to p is unknown
      
      // Search point
      let px = Triangulation.coords[2 * p], py = Triangulation.coords[2 * p + 1]
      
      // Search all triangles
      repeat {
        // Start by getting the orientation of each edge in the current triangle
        let a0 = a - a % 3
        let a1 = a0 + (a + 1) % 3
        let a2 = a0 + (a + 2) % 3
        
        //
        // Is the point p inside the triangle       i
        // This is actually just [e, i, q]         /a\
        //                                        q---e
        //
        
        // Save the vertices
        let q = vertices[a]
        let e = vertices[a1]
        let i = vertices[a2]
        
        // And their co-ordinates
        let ex = Triangulation.coords[2 * e], ey = Triangulation.coords[2 * e + 1]
        let qx = Triangulation.coords[2 * q], qy = Triangulation.coords[2 * q + 1]
        let ix = Triangulation.coords[2 * i], iy = Triangulation.coords[2 * i + 1]
        
        // Orientations?
        let o  = orient(qx, qy, ex, ey, px, py)
        let o1 = orient(ex, ey, ix, iy, px, py)
        let o2 = orient(ix, iy, qx, qy, px, py)
        
        // If none are negative we have found the triangle
        // So let's look for negative orientations
        var found = true
        
        // Edge a1 (q->i)
        if o1 < 0 {
          // Not found
          found = false
          
          // Add the opposing half edge to the stack
          let ha1 = halfEdges[a1]
          if BoundaryEdge < ha1 {
            // Was triangle ha1 already searched?
            let h0 = 3 * (ha1 / 3)
            if !searchedTriangles.contains(h0) { searchStack.append(ha1) }
          }
        } else if 0 == o1 {
          // Is the vertex p between the two endpoints (e, i)
          if isBetween(first: e, second: i, middle: p) {
            onEdge = a1
            return onEdge
          }
          found = false
        }
        
        // Edge a2 (i->e)
        if o2 < 0 {
          // Not found
          found = false
          
          // Add the opposing half edge to the stack
          let ha2 = halfEdges[a2]
          if BoundaryEdge < ha2 {
            // Was triangle ha1 already searched?
            let h0 = 3 * (ha2 / 3)
            if !searchedTriangles.contains(h0) { searchStack.append(ha2) }
          }
        } else if 0 == o2 {
          // Is the vertex p between the two endpoints (i, q)
          if isBetween(first: i, second: q, middle: p) {
            onEdge = a2
            return onEdge
          }
        
          found = false
        }
        
        // Need to consider a
        if o < 0 {
          // Not found
          found = false
          
          // Special case of first call
          if hint == a {
            // Add the opposing half edge to the stack
            // First case so not searched
            let ha = halfEdges[a]
            if BoundaryEdge < ha { searchStack.append(ha) }
          }
        } else if found {
          // If no orientations were negative all done
          // If one orientaton is zero
          if 0 == o {
            // Is the vertex p between the two endpoints (q, e)
            if isBetween(first: q, second: e, middle: p) {
              onEdge = a
              return onEdge
            }
            
            found = false
          }
          return a0
        } // End of findTriangle
        
        // This triangle was searched
        searchedTriangles.insert(a0)
        
        // Get next stack entry
        a = searchStack.popLast() ?? EmptyEdge
        if EmptyEdge == a {
          throw triangulationError.initError("Couldn't find a triangle that encloses vertex \(p) => (\(px), \(py))")}
      } while true
    }
    
    // Use a nested function - this can be unwrapped
    let t = try! findTriangle(edge: hint, vertex: p)
    if EmptyEdge == onEdge {
      return addVertex(inside: t, vertex: p)
    }
    return addVertex(on: onEdge, vertex: p)
  }

 
  //
  mutating func flipEdges(_ e: Int, list flipList: inout Array<Array<Int>>?,
                          flipFour allEdges:Bool = false,
                          repeat flipMany:Bool = true) throws -> Int {
    var a:Int = e

    // Maintains a stack of edges for flipping
    // If an edge needs flipping adds it to the stack - retains edge
    // The order edges are popped off  the stack is important
    // The last popped edge will always end up as being a with
    // edge a2 connecting vertex i->q
    
    // If multiple flips occur a connected region of flipped triangles
    // is created; a perimeter of replaced edges could be recorded
    // So first flip is which creates perimeter [a,b1]
    // If a is flipped into another triangle c this becomes [a, c1, b] and so on
    
    
    // Initially flipList[a] = nil
    var a2:Int
    flipEdge: repeat {
      // A subtlety here is that if an edge [a] is flipped
      // then no new edge is popped from the stack so the edge [a] (the old [b2])
      // is tested again - so in fact two final edges ([a] and [b1]) are checked
      let b:Int = halfEdges[a]

      /* if the pair of triangles doesn't satisfy the Delaunay condition flip it
       *
       *    Input edge is a
       *
       *
       *            n                     n
       *          /| |\                  /  \
       *      a1 / | | \b2            a1/    \a
       *        /  | |  \              /      \
       *       /   | |   \    flip    /___a2___\
       *      i   a| |b  p     =>    i__________p
       *       \   | |   /            \   b2   /
       *        \  | |  /              \      /
       *       a2\ | | /b1             b\    /b1
       *          \| |/                  \  /
       *            q                     q
       *
       *
       */

      // Old vertices are [a] triangle => [a, a1, a2] halfEdges => [q, n, i]
      //               and [b] triangle => [b, b1, b2] halfEdges => [n, q, p]

      // New vertices are [a] triangle => [a, a1, a2] halfEdges => [p, n, i]
      //               and [b] triangle => [b, b1, b2] halfEdges => [i, q, p]
      //
      // So changes    are vertices[a]  => p
      //                   vertices[b]  => i
      //                   halfEdges[a]  <=> halfEdges[b2]
      //                   halfEdges[a2] <=> halfEdges[b2]
      //                   halfEdges[b]  <=> halfEdges[a2]

      // No way to know which of the [0, 1, 2] edges of b were cut in the projection
      // so complex remainder trickery to get the edge ids

      // Get the edges associated with each triangle
      // Triangle a - need a2 for special case
      let a0 = a - a % 3
      a2 = a0 + (a + 2) % 3

      // First case is if (a <-> b) is the convex hell edge
      // In this case no flip can occur
      if b <= BoundaryEdge {
        // If no edges left on stack all edges are locally delaunay
        // Get the next pending edge
        let x = Triangulation.edgeStack.popLast() ?? FinishedEdge
        
        // All done if this is finished
        if x == FinishedEdge { return a2 }
        
        // Next flip
        a = x
        continue flipEdge
      }

      // Finish Triangle a & Triangle b - don't actually need b1 unless added to the stack
      let a1 = a0 + (a + 1) % 3
      let b0 = b - b % 3
      let b2 = b0 + (b + 2) % 3

      // Now the four points
      let n = vertices[a1]
      let i = vertices[a2]
      let q = vertices[a]
      let p = vertices[b2]

      // We need to flip the edges a-b if the point p is inside the circum-circle
      // of the triangle [n, i, q]

      // Use Shewchuk's robust version of inCircumCircle
      let illegal = inCircumCircle(
        ax: Triangulation.coords[2 * n], ay: Triangulation.coords[2 * n + 1],
        bx: Triangulation.coords[2 * i], by: Triangulation.coords[2 * i + 1],
        cx: Triangulation.coords[2 * q], cy: Triangulation.coords[2 * q + 1],
        px: Triangulation.coords[2 * p], py: Triangulation.coords[2 * p + 1])

      // We have to flip this edge
      if illegal {
        vertices[a] = p
        vertices[b] = i

        // The edges a & b are swapped with b2 & a2
        // In general need to see if either a2 or b2 were boundary edges
        // In  the case of one-direction searching (normal) only b2
        let hb2 = halfEdges[b2]
        
        // Recording flips?
        if nil != flipList {
          if hb2 <= BoundaryEdge {
            // This means that a replaces b2
            flipList!.append([a, b2])
          }
 
          // Add latest edge to flip list - a is already added
          if allEdges && (halfEdges[a2] <= BoundaryEdge) {
            flipList!.append([b, a2])
          }
        }
            
        // Link the half edges
        //   halfEdges[a]  <=> halfEdges[b2]
        //   halfEdges[a2] <=> halfEdges[b2]
        //   halfEdges[b]  <=> halfEdges[a2]
        link(a, hb2)
        link(b, halfEdges[a2])
        link(a2, b2)

        // Add the edge [b1] to the stack
        // The edge [a] is not popped so will be tested again
        Triangulation.edgeStack.append(b0 + (b + 1) % 3) // b1
        if allEdges {
          Triangulation.edgeStack.append(b) // b
          Triangulation.edgeStack.append(a1) // a1
        }
      } else {
        // Only do this when no flip happened
        
        // If no edges left on stack all edges are locally delaunay
        let x = Triangulation.edgeStack.popLast() ?? FinishedEdge
        
        // All done if this is finished
        if (x == FinishedEdge) { return a2 }
        
        // Next flip
        a = x
      }
    } while flipMany // End flipEdge
    return a2
  }
  
  // New version
  // Draw the edges of the Voronoi and/or Delaunay structure
  func vtkEdges(label s:String, type flag:String = "Voronoi") -> String {
    let labelPoints = (flag == "Delaunay")

    // Which points are used?
    var usedEdges = 0
    var usedIndex: Dictionary<Int, Int> = [:]
    
    // A bounding box - two virtual points
    let dx = boundingBox[1][0] - boundingBox[0][0]
    let dy = boundingBox[1][1] - boundingBox[0][1]
    
    let xmin = boundingBox[0][0] - 0.05 * dx
    let ymin = boundingBox[0][1] - 0.05 * dy
    let xmax = boundingBox[1][0] + 0.05 * dx
    let ymax = boundingBox[1][1] + 0.05 * dy
    
    //
    var vtkString = "# vtk DataFile Version 2.0\n"
    vtkString += "Edges \(s)\n"
    vtkString += "\nASCII\nDATASET UNSTRUCTURED_GRID\n"
    
    // Which points and edges are in use
    var pString = "", eString = "", cString = ""
    usedVertices: for (e, otherEdge) in halfEdges.enumerated()  {
      if otherEdge != UnusedEdge {
        if emptyTriangles.contains(3 * (e / 3)) { continue usedVertices}
        let value = otherEdge > BoundaryEdge ? 100 : 60
        
        // Is this point used?
        let p = vertices[e]
        
        // Is it a new point?
        if nil == usedIndex[p] {
          usedIndex[p] = usedIndex.count
          pString += String(format: "\n%14.7f %14.7f 0.0", Triangulation.coords[2 * p], Triangulation.coords[2 * p + 1])
          if labelPoints { cString += "\n\(p)" }
        }
        
        // The edges
        let e0 = 3 * (e / 3)
        let f = e0 + (e + 1) % 3

        let q = vertices[f]
        if nil == usedIndex[q] {
          usedIndex[q] = usedIndex.count
          pString += String(format: "\n%14.7f %14.7f 0.0", Triangulation.coords[2 * q], Triangulation.coords[2 * q + 1])
          if labelPoints { cString += "\n\(q)" }
        }
        
        // Most edges occur twice
        if otherEdge <= BoundaryEdge {
          usedEdges += 1
          eString += String(format: "\n2 %d %d", usedIndex[p]!, usedIndex[q]!)
          if !labelPoints { cString += "\n\(value)" }

        } else if p > q {
          usedEdges += 1
          eString += String(format: "\n2 %d %d", usedIndex[p]!, usedIndex[q]!)
          if !labelPoints { cString += "\n\(value)" }
        }
      }
    }
    
    // Save used counts
    let indexCount = usedIndex.count
    
    // Get each triangle
    usedIndex.removeAll(keepingCapacity: true)
    if !labelPoints {
      for t in stride(from: 0, through: numberEdges - 1, by: 3) {
        // Each triangle has three half edges
        // Skip empty triangles
        if emptyTriangles.contains(t) { continue }
        
        // We only need to list each circumcentre once
        let p = triangleCentre(inside: t)
        
        // Save co-ordinates
        pString += String(format: "\n%14.7f %14.7f 0.0", p[0], p[1])
        
        // Don't worry about duplicated points
        usedIndex[t] = indexCount + usedIndex.count // And this has index i
      }
    }
    
    // Voronoi edges join the centre of one circumcircle to another
    vtkString +=  String(format: "\nPOINTS %04d double", usedIndex.count + indexCount + 2)
    vtkString += pString
    
    // A bounding box - two virtual points
    vtkString += String(format: "\n%14.7f %14.7f 0.0", xmin, ymin)
    vtkString += String(format: "\n%14.7f %14.7f 0.0", xmax, ymax)
    
    // Now the edges - only need to draw half the edges
    if !labelPoints {
      for t in usedIndex.keys {
        // Get each neighbouring triangle
        for i in 0...2 {
          let u = 3 * (halfEdges[t + i] / 3)
          if u > t {
            // Draw half of the edges
            // The neighbouring centre
            if nil != usedIndex[u] {
              eString += "\n2 \(usedIndex[t]!) \(usedIndex[u]!)"
              usedEdges += 1
              cString += "\n0"
            }
          }
        }
      }
    }
    vtkString += String(format: "\nCELLS %04d %04d", usedEdges + 2, 3 * usedEdges + 4)
    vtkString += eString
    
    // A bounding box - two virtual cells
    vtkString += String(format: "\n1 %d", usedEdges)
    vtkString += String(format: "\n1 %d", usedEdges + 1)
    
    // the cell types
    vtkString += String(format: "\nCELL_TYPES %04d", usedEdges + 2)
    for _ in 0..<usedEdges {
      vtkString += String(format: "\n3")
    }
    vtkString += String(format: "\n2\n2")
    
    // cell data
    // point data (experimental)
    if labelPoints {
      vtkString += String(format: "\nPOINT_DATA %04d", indexCount + 2)
    } else {
      vtkString += String(format: "\nCELL_DATA %04d", usedEdges + 2)
    }
    
    vtkString += "\nSCALARS edge int 1"
    vtkString += "\nLOOKUP_TABLE default"
    vtkString += cString
    vtkString += "\n75\n75"
    
    return vtkString
  }
  
  // Show the Delaunay triangles
  mutating func vtkDelaunay() -> String {
    // We are done!
    var vtkString = "# vtk DataFile Version 2.0\n"
    vtkString += "Delaunay Mesh \(Zone.name!)\n"
    vtkString += "\nASCII\nDATASET UNSTRUCTURED_GRID\n"
    
    // No editing - just give all the points
    vtkString +=  String(format: "\nPOINTS %04d double", Triangulation.pointCount)
    for p in 0..<Triangulation.pointCount {
      vtkString += String(format: "\n%.7f %.7f 0", Triangulation.coords[2 * p], Triangulation.coords[2 * p + 1])
    }
    
    // These are going to be filled triangles
    let n = (numberEdges / 3) - emptyTriangles.count
    vtkString += String(format: "\nCELLS %04d %04d", n, 4 * n)
    // Get each triangle in turn
    for t in stride(from: 0, through: numberEdges - 1, by: 3) {
      // Each triangle has three half edges
      // Skip empty triangles
      if emptyTriangles.contains(t) { continue }
      
      // List the edges
      vtkString += "\n3 \(vertices[t]) \(vertices[t + 1]) \(vertices[t + 2])"
    }
    
    // the VTK cell types (triangles!)
    vtkString += String(format: "\nCELL_TYPES %04d", n)
    for _ in 0..<n {
      vtkString += String(format: "\n5")
    }
    
    // They are filled so we need some cell data - just the cell id?
    vtkString += String(format: "\nCELL_DATA %04d", n)
    
    vtkString += "\nSCALARS properties int 1"
    vtkString += "\nLOOKUP_TABLE default"
    for i in 0..<n {
      vtkString += "\n\(i)"
    }
    
    return vtkString
  }
  
  // Get a list of Voronoi zones
  mutating func toVoronoi() -> Int {
    // Create triangleCentres
    // and zoneList
    var usedVertices = Set<Int>()
        
    // All Voronoi cells are convex...
    // So simply get all the constituent internal vertices and find their
    // bounding polygons
    var zoneCount = 0
    var e = 0
    eachEdge: while e < halfEdges.count {
      let otherEdge = halfEdges[e]
      if otherEdge > BoundaryEdge {
        // Skip empty Triangles
        let e0 = 3 * (e / 3)
        if emptyTriangles.contains(e0) {
          e += 3; continue eachEdge
        }
        
        // We have a candidate vertex - but does it lie on a loop?
        let v = vertices[e]
        
        // Is this a new vertex?
        if usedVertices.contains(v)  { e += 1; continue eachEdge }
        
        // A new vertex
        usedVertices.insert(v)
        
        // Now get a new polygon
        // Ok we can get each circumcentre
        var polygon = Array<Int>()
        var a = e
        repeat {
          let i  = a / 3
          let a0 = 3 * i
          let a2 = a0 + (a + 2) % 3
          
          // Save index
          // If this point is already computed add it to  the polygon
          var p = triangleCentres[i]
          if nil == p {
            // Compute the centre, save it and link it
            p = triangleCentre(inside: a0)
            triangleCentres[i] = p
            if i >= triangleCount {
              triangleCount = i + 1
            }
          }
          
          // Polygons can have degenerate entries!!
          // The points are mapped to vertices
          checkDuplicate: for j in polygon {
            let q = triangleCentres[j]!
            let x = dist(q[0], q[1], p![0], p![1])
            
            // Is this a matching point?
            if x < squaredThreshold {
              p = nil
              break checkDuplicate
            }
          }
          
          // if not a duplicate
          if nil != p {
            polygon.append(i)
          }
          
          // Get next edge
          a = halfEdges[a2]
        } while a != e && a > BoundaryEdge
        
        // The zone list - suppress external facets
        if a > BoundaryEdge {
          zoneList.append(polygon)
          
          // Each zone has an associated vertex
          zoneVertex.append(v)
          
          // Need to sum the number of entries (One per index plus
          zoneCount += polygon.count + 1
        }
      }
      
      // Next edge
      e += 1
    }
        
    // A modified zone count used by VTK
    return zoneCount
  }
 
  // Write a VTK formatted version (showing zone cells!)
  mutating func vtkVoronoi(_ zoneCount:Int) -> String {
    // We are done!
    var vtkString = "# vtk DataFile Version 2.0\n"
    vtkString += "Voronoi Diagram \(Zone.name!)\n"
    vtkString += "\nASCII\nDATASET UNSTRUCTURED_GRID\n"
    
    // The points are the centres of the circumcircles
    vtkString +=  String(format: "\nPOINTS %04d double", triangleCount)
    
    // Average point for empty centres
    let centreString = String(format: "\n%.1f %.1f 0", 0.5 * (boundingBox[0][0] + boundingBox[1][0]), 0.5 * (boundingBox[0][1] + boundingBox[1][1]))
    
    for i in 0..<triangleCount {
      let p = triangleCentres[i]
      if nil != p {
        vtkString += String(format: "\n%.7f %.7f 0", p![0], p![1])
      } else {
        vtkString += centreString
      }
    }
    
    // These are going to be filled cells - so we need a complete polygon
    vtkString += String(format: "\nCELLS %04d %04d", zoneList.count, zoneCount)
    for polygon in zoneList {
      vtkString += "\n\(polygon.count)"
      for p in polygon {
        vtkString += " \(p)"
      }
    }
    
    // the cell types
    vtkString += String(format: "\nCELL_TYPES %04d", zoneList.count)
    for _ in zoneList {
      vtkString += String(format: "\n7")
    }
    
    // cell data
    vtkString += String(format: "\nCELL_DATA %04d", zoneList.count)
    
    vtkString += "\nSCALARS properties double 1"
    vtkString += "\nLOOKUP_TABLE default"
    for (_, v) in zoneVertex.enumerated() {
      let z = Zone.propertyList[Triangulation.properties[v]].tag ?? Zone.propertyList[Triangulation.properties[v]].numberPoints
      vtkString += "\n\(z ?? 0)"
    }
    
    return vtkString
  }
  
  // Create a yaml formatted zone file
  func toZone() -> String {    
    // We are done!!
    var string = "# Yaml File\ntype: Vorobox\n"

      string += "control:\n"
      string += "  showMe: \(showMe)\n"
      string += "  iteration: \(Zone.control.iteration! + 1)\n"

    // Geometry
    string += "objects:\n"
    string += "  \(Zone.name!):\n"
    string += "    - type: GeometryCollection\n"
    string += "      geometries:\n"
    string += "        - multipolygon: ["
    for i in 0..<zoneList.count {
      string += "[[\(i)]], "
    }
    string += "]\n"
    
    // Now the arcs
    string += "arcs: [\n"
    for (i, polygon) in zoneList.enumerated() {
      string += "[&_\(i) "
      for p in polygon {
        let c = triangleCentres[p]!

        string += "[\(c[0]), \(c[1])],"
      }
      string += " *_\(i)],\n"
    }
    string += "]\n"
    
    return string
  }
  

  // circumcentre of a triangle
  func triangleCentre(inside e:Int) -> [Double]  {
    return circumCentre(Triangulation.coords[2 * vertices[e]], Triangulation.coords[2 * vertices[e] + 1],
                        Triangulation.coords[2 * vertices[e + 1]], Triangulation.coords[2 * vertices[e + 1] + 1],
                        Triangulation.coords[2 * vertices[e + 2]], Triangulation.coords[2 * vertices[e + 2] + 1])
  }

  // Logging
  internal func showVoronoi(label s:String, type flag: String = "Voronoi", force show: Bool = false) {
    showCount += 1
    let i = showCount % showLimit
    
    // Debugging - More information
    if show || showMe > 2  {
      // And a vtk file - voronoi
      let folder = try! Folder(path: OutputFolder)
      let filename = "edges_" + s + String(format:"%04d.vtk", i)
      let file = try! folder.createFile(named: filename)
      try! file.write(vtkEdges(label: "\(i)", type: flag))
    }
  }
}
  

/* Helper Functions
 dist:
 inCircumCircle:
 circumCentre:
 circumRadius:

 */
func dist(_ ax: (Double), _ ay: (Double),
          _ bx: (Double), _ by: (Double)) -> Double {
  let dx = ax - bx
  let dy = ay - by
  return dx * dx + dy * dy
}

// Convenience method - use indices as a wrapper
func insideCircumCircle(vertex p:Int, _ n:Int, _ i:Int, _ q:Int) -> Bool {  
  // Use Shewchuk's robust version of inCircumCircle
  return inCircumCircle(ax: Triangulation.coords[2 * n], ay: Triangulation.coords[2 * n + 1],
                        bx: Triangulation.coords[2 * i], by: Triangulation.coords[2 * i + 1],
                        cx: Triangulation.coords[2 * q], cy: Triangulation.coords[2 * q + 1],
                        px: Triangulation.coords[2 * p], py: Triangulation.coords[2 * p + 1])
}

// Is point p inside the triangle abc
// Documentation - Shewchuk
func inCircumCircle(ax: (Double), ay: (Double),
              bx: (Double), by: (Double),
              cx: (Double), cy: (Double),
              px: (Double), py: (Double)) -> Bool {
  let adx = ax - px
  let ady = ay - py
  let bdx = bx - px
  let bdy = by - py
  let cdx = cx - px
  let cdy = cy - py

  let bdxcdy = bdx * cdy
  let cdxbdy = cdx * bdy
  let alift = adx * adx + ady * ady

  let cdxady = cdx * ady
  let adxcdy = adx * cdy
  let blift = bdx * bdx + bdy * bdy

  let adxbdy = adx * bdy
  let bdxady = bdx * ady
  let clift = cdx * cdx + cdy * cdy

  let det = alift * (bdxcdy - cdxbdy)
    + blift * (cdxady - adxcdy)
    + clift * (adxbdy - bdxady)

  let permanent = (abs(bdxcdy) + abs(cdxbdy)) * alift
    + (abs(cdxady) + abs(adxcdy)) * blift
    + (abs(adxbdy) + abs(bdxady)) * clift
  let errbound = circleThreshold * permanent
  
  // Positive - anticlockwise
  // Negative - clockwise
  return det > errbound
}

// Get the circum-radius squared
func circumRadius(_ px: Double, _ py: Double,
                  _ qx: Double, _ qy: Double,
                  _ rx: Double, _ ry: Double) -> Double {
  // Denominator
  let c = circumCentre(px, py, qx, qy, rx, ry)
  let dx = c[0] - rx, dy = c[1] - ry
  return dx * dx + dy * dy // r2
}

// Get the centre of the circumcircle
func circumCentre(_ px: (Double), _ py: (Double),
                  _ qx: (Double), _ qy: (Double),
                  _ rx: (Double), _ ry: (Double)) -> [Double] {
  // Check orientation
  let o2d = orient(px, py, qx, qy, rx, ry)
  if o2d == 0 { return [Double.infinity, Double.infinity]}
  
  let pxrx = px - rx, pyry = py - ry
  let qxrx = qx - rx, qyry = qy - ry
  let prdx = pxrx * pxrx + pyry * pyry
  let qrdx = qxrx * qxrx + qyry * qyry
  let detx = qyry * prdx - pyry * qrdx
  let dety = pxrx * qrdx - qxrx * prdx
  
  return [rx + 0.5 * detx / o2d, ry + 0.5 * dety / o2d]
}

// Get the centre of the incircle
func inCentre(_ ax: (Double), _ ay: (Double),
              _ bx: (Double), _ by: (Double),
              _ cx: (Double), _ cy: (Double)) -> [Double] {
  // For efficiency keep as a separate code bock
  // Lengths of edges opposite the vertices c, a, b
  let c = dist(ax, ay, bx, by).squareRoot()
  let a = dist(bx, by, cx, cy).squareRoot()
  let b = dist(cx, cy, ax, ay).squareRoot()
  let abc = a + b + c // The perimeter
  return [(a * ax + b * bx + c * cx) / abc, (a * ay + b * by + c * cy) / abc]
}

// return 2d orientation sign if we're confident in it through J. Shewchuk's error bound check
//
// Is a point c to the left or right of a vector a->b
// Or more symmetrically are the points (a->b->c) in an anti-clockwise order? +ve yes -ve no
// The magnitude of this quantity is twice the triangle area
//
/**
 [Robust Predicates](http://www.cs.cmu.edu/~quake/robust.html)
 */
func orient(_ ax: (Double), _ ay: (Double),
                  _ bx: (Double), _ by: (Double),
                  _ cx: (Double), _ cy: (Double)) -> Double {

  // Points a, b, c
  // Form a triangle

  // The left and right determinants
  let detRight = (ay - cy) * (bx - cx)
  let detLeft  = (ax - cx) * (by - cy)

  // We only need to determine the sign +ve/zero/-ve
  // Positive - anticlockwise
  // Zero     - colinear
  // Negative - clockwise
  let det = detLeft - detRight

  return abs(det) > orientThreshold * abs(detLeft + detRight) ? det : 0
}

// A wrapper
//func orient(_ ax: (Double), _ ay: (Double),
//            _ bx: (Double), _ by: (Double),
//            _ cx: (Double), _ cy: (Double)) -> Bool {
//  // Just a wrapper
//  return orient(ax, ay, bx, by, cx, cy) > 0
//}

func isBetween(first i:Int, second j:Int, middle k:Int) -> Bool {
  // Is the point k between the points i & j
  let ax = Triangulation.coords[2 * i]
  let bx = Triangulation.coords[2 * j]
  
  // Try abscissae first
  if ax != bx {
    let cx = Triangulation.coords[2 * k]
    return (ax - cx) * (bx - cx) <= 0
  }
  
  // Use the ordinates
  let ay = Triangulation.coords[2 * i + 1]
  let by = Triangulation.coords[2 * j + 1]
  let cy = Triangulation.coords[2 * k + 1]
  return (ay - cy) * (by - cy) <= 0
}

// Reflect the point c in the line from a to b
func reflect(_ cx:Double, _ cy:Double, _ ax:Double, _ ay:Double, _ bx:Double, _ by:Double) -> Int {
  // It's useful to compute some properties of this mirror line
  let count = Triangulation.coords.count / 2
  let dx  = bx - ax
  let dy  = by - ay
  let dx2 = dx * dx + dy * dy // dist(....)
  let d1  = (dx * dx - dy * dy) / dx2
  let d2  = 2 * dx * dy / dx2
  
  // The reflected point co-ordinates
  let rx  = d1 * (cx - ax) + d2 * (cy - ay) + ax
  let ry  = d2 * (cx - ax) - d1 * (cy - ay) + ay
  Triangulation.coords.append(rx)
  Triangulation.coords.append(ry)
  Triangulation.pointCount += 1
  
  // Vertex id
  return count
}

// Obtain a point a unit length away from c along a vector perpendicular to the vector
// from a to b after an anticlockwise rotation
func perpendicular(scale r2:Double, _ cx:Double, _ cy:Double, _ ax:Double, _ ay:Double, _ bx:Double, _ by:Double) -> Array<Double> {
  // It's useful to compute some properties of this mirror line
  let dx  = bx - ax
  let dy  = by - ay
  let f = (r2 / (dx * dx + dy * dy)).squareRoot() // dist(....)
    
  // The rotated point co-ordinates
  let rx  = cx - f * dy
  let ry  = cy + f * dx
  
  // The point
  return [rx, ry]
}

/* Deprecated functions to be replaced in final clean up */

// Near zero
func isZero(_ x: Double) -> Bool {
  return (x <= Epsilon) && (x >= -Epsilon)
}
