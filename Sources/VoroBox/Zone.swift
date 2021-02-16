//  VoroBox
//
//  Zone.swift
//
//  Created by Rob Whitehurst on 26/5/20.
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
//
import Foundation

// Properties
struct Properties: Codable {
  // Zone properties
  var density:Double = 1
  var name:String?
}

// Control 
struct Control: Codable {
  // Zone properties
  var id:Int = FirstZone
  var iteration:Int = 0
  var distinct:Bool = true
  var showMe:Int = 0
}

extension Control {
  init(from c:StoredZones.Control) {
    id = c.id ?? FirstZone
    distinct = c.distinct ?? true
    showMe = c.showMe ?? 0
    iteration = c.iteration ?? 0
  }
}

enum zoneError: Error {
  case initError(String)
}

// Debugging
internal var showMe = 0

// The Zone structure
struct Zone: Codable {
  // A global name
  static var name:String? = nil
  
  // This should be a Triangulation Boolean
  static var hullConforming: Bool?
      
  // Arc vertices (when using topojson style zone files)
  static var arcVertices = Array<Array<Int>>()

  // Zone properties 
  static var propertyList = Array<Properties>()
  
  // Global properties
  static var control:Control = Control()
  static var origin:Array<Double> = [0, 0]
  static var scale:Array<Double> = [1, 1]
  
  // Instance properties
  var origin:Array<Double> = [0, 0]
  var scale:Array<Double> = [1, 1]

  // The polygon bounding the zone
  // along with its size
  var polygon = Array<String>()
  var polygonCount = 0

  // Convex zones
  var convexZones = Array<Zone>()
  
  // Indices are associated with each instance
  var zoneIndices = Array<Int>()
  
  // Each zone can have zero or more holes
  var holeIndices = Array<Array<Int>>()
    
  // A code which identifies whether a zone is delaunay or voronoi conforming
  // This is deprecated
  var code: Int = 1

  // Index to find the properties associated with this zone
  var propertyIndex: Int = 0
}

// The actual zone structure ...
extension Zone {
  // Read in a zone defined in a data file
  init(usingData stored: StoredZones.Zone, with arcs:Array<Array<Array<Double>>>) throws {
    var inputZones = Array<Zone>()
    
    // Global controls 
    var distinct = Zone.hullConforming ?? true
    showMe = Zone.control.showMe 
    
    // Set the zone code
    code = 2 * Zone.control.id

    // Update the code depending on if this zone is distinct or not
    if distinct { code += 1 } // Code is odd for conforming edges
    
    // Set code - when even the zone is Voronoi conforming (the default)
    Zone.control.id += 1
    
    // Make space for the arcs
    Zone.arcVertices = Array(repeating: [], count: arcs.count)
    
    // Initially we are reading the stored data
    // // Link properties to zones 
    // if let p = stored.properties {
    //   Zone.propertyList.append(Properties(from: p))
    // } else if let list = stored.multiproperties {
    //   for p in list {
    //    Zone.propertyList.append(Properties(from: p))
    //   }
    // }
        
    // We are given the co-ordinates
    // Only  the first array is needed to construct the boundary
    //
    // This might be a multi-part zone
    // Properties are always noted for every zone 
    var propertyIndex = Zone.propertyList.count
    if let polygons = stored.polygon {
      // Make polygons arrays of vertices
      // Each polygon is a [[#boundary_arcs#],[#hole-arcs#],...]
      
      // // Define a zone using a polygon
      // if let p  = stored.properties {
      //   Zone.propertyList.append(Properties(from: p))
      // } else {
      //   // Empty properties
      //   Zone.propertyList.append(Properties())
      // }
      inputZones.append(Zone(copy: self, using: [polygons], index: 0, and: arcs, data: stored))
    } else if let multipolygons = stored.multipolygon {
      // // Link the properties 
      // if let list = stored.multiproperties {
      //   if list.count != multipolygons.count {
      //     throw triangulationError.initError("MultiProperties list length \(list.count) does not match the multipolygon list length \(multipolygons.count)")
      //   }

      //   for p in list {
      //    Zone.propertyList.append(Properties(from: p))
      //   }
      // } else {
      //    Zone.propertyList.append(Properties(from: p))
      // }

      // This is quite straightforward now
      for i in 0..<multipolygons.count {
        inputZones.append(Zone(copy: self, using: multipolygons, index: i, and: arcs, data: stored))
      }


    }
    
    // Save the properties 


    // Explicitly specified points
    var specifiedPoints = Array<Array<Double>>()
    if let point = stored.point {
      specifiedPoints.append(point)
    } else if let points = stored.multipoint {
      specifiedPoints.append(contentsOf: points)
    }
    
    if showMe > 0 {
      print("Processing \(inputZones.count) zones")
    }
    for i in 0..<inputZones.count {
      var z = inputZones[i]
      if !z.zoneIndices.isEmpty {
        // Now need to make a list of convex zones
        try! z.clipZone()
      } // Clipped the zone into convex polygons
    
      // Convex zones include z
      z.convexZones.append(z)
      
      // Add specified points to the zone
      for i in 0..<z.convexZones.count {
        // Add the points - once one is added it need not be considered again
        specifiedPoints = z.convexZones[i].addPointsToZone(specifiedPoints)
      }
      
      // The final convex zone is simply z
      convexZones.append(contentsOf: z.convexZones)
      //convexZones.append(z)
      
      if showMe > 0 {
        print("Added \(convexZones.count) convex zones")
      }
    }
  } // init
 
  mutating func addPointsToZone(_ pointList: Array<Array<Double>>) -> Array<Array<Double>> {
    // Now we can allocate the points (if any) that are associated with each convex zone
    var newList = Array<Array<Double>>()
    
    var minX = boundingBox[0][0]
    var minY = boundingBox[0][1]
    var maxX = boundingBox[1][0]
    var maxY = boundingBox[1][1]
    
    // Add defined points
    nextPoint: for p in pointList {
      // The point
      let x = p[0], y = p[1]
      
      // Is this inside this zone?
      if polygonCount == 0 || isInside(boundary:zoneIndices, size: polygonCount, query_x: x, query_y: y) {
        // Exclude points that are inside any hole
        for h in holeIndices {
          if isInside(boundary: h.reversed(), size: h.count, query_x: x, query_y: y) {
            // Inside a hole
            continue nextPoint
          }
        }
        
        // Save this point in the coordinate list
        // Add it to the coordinate list
        Triangulation.coords.append(x)
        Triangulation.coords.append(y)
        
        // These are given a zone code
        Triangulation.code.append(code)
        
        // Yes this point is inside this zone
        zoneIndices.append(Triangulation.pointCount)
        
        // One extra point
        Triangulation.pointCount += 1
        
        // Compute bounds - make no assumption about labelled points
        if (x < minX)  {minX = x}
        if (y < minY)  {minY = y}
        if (x > maxX)  {maxX = x}
        if (y > maxY)  {maxY = y}
        
        // Can only be in one zone
      } else {
        // Check p in another list
        newList.append(p)
      }
    }
    
    // All done if no perimeter or bounding box
    if polygonCount < 2 { return newList}
    
    // Add points defined by a perimeter
    var i = 1
    while i < polygonCount {
      let x = Triangulation.coords[2 * zoneIndices[i]], y = Triangulation.coords[2 * zoneIndices[i] + 1]
      
      // Compute bounds
      if (x < minX)  {minX = x}
      if (y < minY)  {minY = y}
      if (x > maxX)  {maxX = x}
      if (y > maxY)  {maxY = y}
      
      // next point
      i += 1
    }
    
    // Update the overall bounding box
    if (boundingBox[0][0] > minX) { boundingBox[0][0] = minX }
    if (boundingBox[0][1] > minY) { boundingBox[0][1] = minY }
    if (boundingBox[1][0] < maxX) { boundingBox[1][0] = maxX }
    if (boundingBox[1][1] < maxY) { boundingBox[1][1] = maxY }
    
    // Now add random points
    // Update zone count
    Triangulation.pointCount = Triangulation.coords.count / 2

    // This is a clumsy make do
    // Can be replaced once other components working
    let density:Double = Double.random(in: 0...100)
    
    // Number of points is density times box area
    var r:Double = density * (maxX - minX) * (maxY - minY)
    r.round(.up)
      
    // Add the points
    var ii = Int(r)
    nextRandomPoint: while ii > 0 {
      // Add points inside  the zone
      let x = Double.random(in: minX...maxX), y = Double.random(in: minY...maxY)
      
      // Append the point
      if isInside(boundary: zoneIndices, size: polygonCount, query_x: x, query_y: y) {
        // Skip any inside a hole
        for h in holeIndices {
          if isInside(boundary: h.reversed(), size: h.count, query_x: x, query_y: y) {
            // Inside a hole
            continue nextRandomPoint
          }
        }
        
        // New point - set coordinates
        Triangulation.coords.append(x)
        Triangulation.coords.append(y)
        
        // This needs the zone code
        Triangulation.code.append(code)
        
        // Set zone index
        zoneIndices.append(Triangulation.pointCount)
        Triangulation.pointCount += 1

        ii -= 1
      }
    } // Add random points
    
    return newList
  }
  
  
  // This copies a zone but replaces the boundary with the specified arcs
  init(copy zone: Zone, using multipolygons:Array<Array<Array<Int>>>, index polyIndex:Int,
       and arcs:Array<Array<Array<Double>>>,
       data stored: StoredZones.Zone) {
    func polyFrom(arcs arcList:Array<Int>, check beforeIndex:Int = 0) -> Array<Int> {
      // Passed in an array listing the arcs to be used
      // This is used to (i) construct the polygon of vertices
      // and (ii) - as a side effect - store the vertices for future use
      // The lists of vertices that defines the polygon,
      // first list is the
      var polyVertices = Array<Int>()

      // Some arcs (boundary arcs in multipolygons) also need
      // to check other earlier boundary arcs in the multipolygon
      // in order to establish if there are duplicated points
      
      // Get each arc in turn
      for a in arcList {
        // This is a list of arc indices
        // Each arc is indexed by a
        // Is this arc read in reverse order?
        let reverse = a < 0 ? true : false
        
        // Get the arc index (ones complement for negative numbers)
        let arcIndex = reverse ? ~a : a
                          
        // If this arc was examined before the vertices will be listed already
        if Zone.arcVertices[arcIndex].isEmpty {
          // The arc is itself a list of points
          // The last point in each arc is the same as the first in the next arc
          // Generate the vertices
          let inputArc = arcs[arcIndex] // An Array<[Double,Double]>
          
          // The points are mapped to vertices
          roundList: for k in 0..<inputArc.count {
            var p = inputArc[k]
            
            // Can reuse vertices - this can only happen in the same arc
            if k > 0 {
              // Is this a matching point?
              // A duplicated point must be in the same arc
              // Get each earlier point in this arc
              if let vertexIndex = inputArc[0...(k-1)].firstIndex(of: p) {
                Zone.arcVertices[arcIndex].append(Zone.arcVertices[arcIndex][vertexIndex])
                continue roundList
              }
            }

            // Do we need to check other arcs?
            if beforeIndex > 0 {
              // Yes 
              for i in 0..<beforeIndex {
                // The polygon boundaries to check 
                let boundaryArcs = multipolygons[i].first!

                // Does p exist on these arcs?
                for b in boundaryArcs {
                  // This is a list of arc indices
                  // Each arc is indexed by b                  
                  // Get the arc index (ones complement for negative numbers)
                  let boundaryIndex = b < 0 ? ~b : b
                  
                  // Is the point in the boundary?
                  if let vertexIndex = arcs[boundaryIndex].firstIndex(of: p) {
                    Zone.arcVertices[arcIndex].append(Zone.arcVertices[boundaryIndex][vertexIndex])
                    continue roundList
                  }
                }
              }
            }
            
            // A new point if none found -
            // Transform the point p
            for j in 0...1 {
              p[j] = origin[j] + scale[j] * p[j]
            }
            
            // Save this point
            Triangulation.coords.append(contentsOf: p)
            
            // Use the zone properties
            Triangulation.code.append(NoZoneCode)
            
            // Save the vertex
            Zone.arcVertices[arcIndex].append(Triangulation.pointCount)
            Triangulation.pointCount += 1
          } // End of each list
        }
          
        // The vertices are now available
        // In building the polygon loops
        // note that the last & first vertices overlap
        if reverse {
          polyVertices.append(contentsOf: Zone.arcVertices[arcIndex].reversed())
        } else {
          polyVertices.append(contentsOf: Zone.arcVertices[arcIndex])
        }
        
        // Drop the last vertex
        polyVertices.removeLast()
      } // end of the arc list
      
      // Return the list of vertices
      return polyVertices
    }
    
    // Add a zone defined by a polygon
    code = zone.code

    // This adds the pointer to one more property 
    propertyIndex = Zone.propertyList.count + polyIndex
  
    // Zone scale (mapped by global scale if reqired)
    // Zones can be referred to a specific origin
    if let o = stored.origin {
      origin  = [o[0] * Zone.scale[0], o[1] * Zone.scale[1]]
    } else {
      origin  = Zone.origin
    }
    
    // And a scale factor can be applied
    if let s = stored.scale {
      scale  = [s[0] * Zone.scale[0], s[1] * Zone.scale[1]]
    } else {
      scale  = Zone.scale
    }

    //
    // Use arcs to define points
    let polygons = multipolygons[polyIndex]

    // For boundary arcs need to consider all the
    // earlier boundary arcs in checking for
    // duplicated vertices
    var polygon = polyFrom(arcs: polygons.first!, check: polyIndex)
     
    // Check ordering - boundaries should be anticlockwise
    if isClockwise(polygon) { polygon.reverse() }

    // The zone
    zoneIndices.append(contentsOf: polygon)

    // Read the holes
    for h in 0..<polygons.count {
      if h > 0 {
        var hole = polyFrom(arcs: polygons[h])

        // Check ordering -  holes should be clockwise
        if !isClockwise(hole) { hole.reverse() }

        // Save this
        holeIndices.append(hole)
      } // End of each hole
    }
    
    // All done
    polygonCount = zoneIndices.count
  }
  
  // Some specialized inits
  init(clip zone: Zone, ear i:Int, order n:Int) {
    // This handles the case of convex zones
    // This generates a triangular zone
    polygonCount = n // Its an n-gon - can be zero
    
    // clipped zones inherit parent code
    code = zone.code
    
    // Needed for initialization
    propertyIndex = zone.propertyIndex
    
    // Copy the perimeter
    // This always relates to the parent instance
    for k in 0..<polygonCount {
      // populate an array of point indices
      // This is the perimeter
      let index = (i + k + zone.polygonCount) % zone.polygonCount
      zoneIndices.append(zone.zoneIndices[index])
    }
  }
  
  // Need to initialize a zone using just a polygon
  init(given perimeter: Array<Int>, code boundaryType:Int = -2) {
    // Lightweight zone
    
    // We have the indices already
    zoneIndices = perimeter
    polygonCount = zoneIndices.count
    code = -boundaryType
    
    // Now need to make a list of convex zones
    do {
      try clipZone()
    } catch {
      fatalError("Couldn't triangulate \(self)\n\(error)")
    }
    
    // The final convex zone is simply self
    convexZones.append(self)
  }
  
//  mutating func addPointsToZone(_ pointList: Array<Int>) -> Array<Int> {
//    // Now we can allocate the points (if any) that are associated with each convex zone
//    var newList = Array<Int>()
//
//    // Add defined points
//    nextPoint: for i in pointList {
//      // Do we have an unattached specified point?
//      let x = Triangulation.coords[2 * i], y = Triangulation.coords[2 * i + 1]
//
//      // Is this inside this zone?
//      if polygonCount == 0 || isInside(boundary:zoneIndices, size: polygonCount, query_x: x, query_y: y) {
//        // Exclude points that are inside any hole
//        for h in holeIndices {
//          if isInside(boundary: h.reversed(), size: h.count, query_x: x, query_y: y) {
//            // Inside a hole
//            continue nextPoint
//          }
//        }
//
//        // Yes this point is inside this zone
//        zoneIndices.append(i)
//
//        // Can only be in one zone
//      } else {
//        newList.append(i)
//      }
//    }
//
//    // All done if no perimeter or bounding box
//    if polygonCount < 2 { return newList}
//
//    // Add points defined by a perimter
//    var i = 1
//    var minX = Triangulation.coords[2 * zoneIndices[0]]
//    var minY = Triangulation.coords[2 * zoneIndices[0] + 1]
//    var maxX = minX
//    var maxY = minY
//
//    while i < polygonCount {
//      let x = Triangulation.coords[2 * zoneIndices[i]], y = Triangulation.coords[2 * zoneIndices[i] + 1]
//
//      // Compute bounds
//      if (x < minX)  {minX = x}
//      if (y < minY)  {minY = y}
//      if (x > maxX)  {maxX = x}
//      if (y > maxY)  {maxY = y}
//
//      // next point
//      i += 1
//    }
//
//    // Now add random points
//    // Update zone count
//    Triangulation.pointCount = Triangulation.coords.count / 2
//    if let density:Double = properties!.randomDensity {
//      // Number of points is density times box area
//      var r:Double = density * (maxX - minX) * (maxY - minY)
//      r.round(.up)
//
//      // Add the points
//      nextRandomPoint: for _ in 0..<Int(r) {
//        // Add points inside  the zone
//        let x = Double.random(in: minX...maxX), y = Double.random(in: minY...maxY)
//
//        // Append the point
//        if isInside(boundary: zoneIndices, size: polygonCount, query_x: x, query_y: y) {
//          // Skip any inside a hole
//          for h in holeIndices {
//            if isInside(boundary: h.reversed(), size: h.count, query_x: x, query_y: y) {
//              // Inside a hole
//              continue nextRandomPoint
//            }
//          }
//
//          // New point - set coordinates
//          Triangulation.coords.append(x)
//          Triangulation.coords.append(y)
//
//          // This needs the zone code
//          Triangulation.code.append(code)
//
//          // Set zone index
//          zoneIndices.append(Triangulation.pointCount)
//          Triangulation.pointCount += 1
//        }
//      }
//    } // Add random points
//
//    return newList
//  }
  
  // Extra geometric primitives for zones
  // Clip the zone
  // Consider holes
  mutating func clipZone() throws {
    var conflictingLoop = Array<Int>()

    // Outer loop manages the holes
    untilEmpty: repeat {
      var conflictsWith = EmptyVertex
      
      // Nested function to find conflicting indices
      func conflicts(with i:Int) -> Array<Int>? {
        //
        //           i----- *
        //          /     .
        //         /   .  (X) Is point (X) inside the triangle [i-1, i, i+1]
        //        / .
        //       *
        //
        // Examine each vertex on each hole
        // to see if it conflicts with this triangle, and if so
        // record the one that approaches to i closest along
        // a line perpendicular to the edge (i+1 -> i-1)
        var listIndices = Array<Int>()
        
        // The vertices associated with triangle;
        let a = zoneIndices[(i + polygonCount - 1) % polygonCount]
        let b = zoneIndices[i]
        let c = zoneIndices[(i + polygonCount + 1) % polygonCount]

        // And their co-ordinates [abc]
        let ax = Triangulation.coords[2 * a], ay = Triangulation.coords[2 * a + 1]
        let bx = Triangulation.coords[2 * b], by = Triangulation.coords[2 * b + 1]
        let cx = Triangulation.coords[2 * c], cy = Triangulation.coords[2 * c + 1]

        // Save double the area
        var max2D = -Double.infinity
        var closest:(Int, Int)? = nil
        
        // Get each hole
        for (h, hole) in holeIndices.enumerated() {
          // Get each point in each hole
          for (index, p) in hole.enumerated() {
            // This is vertex p
            // Is this point inside the triangle [abc]?
            let px = Triangulation.coords[2 * p], py = Triangulation.coords[2 * p + 1]

            // Orientations?
            if orient(ax, ay, bx, by, px, py) >= 0 && orient(bx, by, cx, cy, px, py) >= 0 {
              let o2d = orient(cx, cy, ax, ay, px, py)
              if o2d >= max2D {
                max2D = o2d
                closest = (h, index)
              }
            }
          }
        } // End of scanning each list of vertices on holes
        
        // Did any conflict?
        if max2D >= 0 {
          // Yes so update the list
          // Rotate the list so that p is the first index
          let (h, p) = closest!
          if p > 0 {
            let count = holeIndices[h].count
            listIndices.append(contentsOf: holeIndices[h].dropFirst(p))
            listIndices.append(contentsOf: holeIndices[h].dropLast(count - p))
          } else {
            listIndices = holeIndices[h]
          }
          
          // Remove holeIndices
          holeIndices[h].removeAll()
          holeIndices.remove(at: h)

          // All done
          return listIndices
        }
        
        // No conflicts
        return nil
      }
      
      // Anything less than quadrilateral is convex
      if polygonCount < 4 && holeIndices.isEmpty { return }
      
      // get the initial ear list
      var ears:Array<Int> = try listEars()
      
      // Determine which ears should be clipped
      // Finished when all are ears
      // Keep track of the last clipped zone
      var convexZone:Zone? = nil
      
      // This can be improved
      untilConvex: while polygonCount > 3 && (!holeIndices.isEmpty || ears.count < polygonCount) {

        // For testing turn off reversible clip
        var reversibleClip = true
         
        // Get the ear at the end of the list
        let i = ears.last!
        
        // Clip corner i
        // If no hole indices conflict with the ear go ahead and clip it
        if let c = conflicts(with: i) {
          // Save the conflicting loop and the relevant boundary ear
          conflictingLoop = c
          conflictsWith = i
          break untilConvex
        } else {
          conflictingLoop = []
        }
        
        // New triangular zone (i-1, i, i+1)
        let z = Zone(clip: self, ear: i - 1, order: 3)
        
        // Update this polygon
        zoneIndices.remove(at: i)
        polygonCount -= 1
        
        // Update ear status - but not needed for
        // triangular zones
        if polygonCount > 3 {
          // update ear status - update polygon - n.b. all done if triangular
          // Remove last ear
          _ = ears.popLast()
          
          // Status of points (i - 1) & (new) point i has changed
          // Note that (i - 1) >= 0 since i >= 1
          // Also record if this is reversible clip
          if isDiagonal(first: i - 2, second: i) {
            // Point i - 1 is now an ear
            if ears.last != (i - 1) {
              // Created a new ear
              reversibleClip = false
              ears.append(i - 1)
            }
          } else { // (i - 1) is not an ear
            // Remove ear is this point present
            reversibleClip = false
            if ears.last == (i - 1) {
              _ = ears.popLast()
            }
          }
          
          // Check if (new) point i is an ear
          if isDiagonal(first: i - 1,
                        second: i + 1) {
            // Yes so is the next recorded ear at point i
            // But this might be i == polygonCount
            //
            if i % polygonCount > 0 {
              // No so can simply append the new ear to the list
              ears.append(i)
            } else if ears.first != 0 {
              // If the zeroth point not already an ear
              // then insert an extra ear at the start
              ears.insert(0, at:0)
            }
            
            // Not reversible
            reversibleClip = false
            
          } else { // point [i] is not an ear
            // Remove ear if this point is present
            // But this only happens if [i] is [0]
            // and if ears[0] == 0
            if ears.first == (i % polygonCount) {
              reversibleClip = false
              ears.removeFirst()
            }
          }
        }
        
        // Two ears theorem 
        if ears.count < 2 {
          print("clipZone: For polygon order \(polygonCount) only \(ears.count) - should be at least two ears")
          throw triangulationError.initError("Less than two ears in clipZone - malformed input?")
        }
        // Update convexZone with the clip
        convexZone = nil == convexZone ? z : try! convexZone!.addTriangularZone(triangle: z)
        
        // If a clip is not reversible we can clear the saved zone
        if !reversibleClip {
          // Not reversible so any saved zone is stored
          convexZones.append(convexZone!)
          
          // Update the saved zone
          convexZone = nil
        }
      } // Until convex
      
      // Got to here - still need to save any unattached zone
      if convexZone != nil {
        convexZones.append(convexZone!)
      }
      
      // Any conflicting indices?
      if !conflictingLoop.isEmpty {
        // Need to splice this loop in
        polygonCount += conflictingLoop.count + 2
        
        // Create a link
        zoneIndices.insert(contentsOf: [zoneIndices[conflictsWith], conflictingLoop.first!], at: conflictsWith)
        
        // Add the loop
        if conflictsWith > 0 {
          zoneIndices.insert(contentsOf: conflictingLoop, at: conflictsWith + 1)
        } else {
          zoneIndices.append(contentsOf: conflictingLoop)
        }
      }

    } while !conflictingLoop.isEmpty

  } // End clips
  
  mutating func addTriangularZone(triangle t:Zone) throws -> Zone {
    // Given a trangular zone t add it to
    // another convex zone self
    // The connecting vertices
    let p = t.zoneIndices[1], q = t.zoneIndices[2]
    
    // Find p in self
    if let i = zoneIndices.firstIndex(of: p) {
      if zoneIndices[(i + polygonCount - 1) % polygonCount] == q {
        // One extra entry on the perimeter
        zoneIndices.insert(t.zoneIndices[0], at: i)
        
        // OK we just grew the polygon
        polygonCount += 1
        
        return self
      }
    }
    
    // Oops
    throw zoneError.initError("Zone \(Zone.name!) can't find matching vertices to join two zones")
  }
  
  //
  // listEars :  Compute the ear status of all points
  func listEars() throws -> Array<Int> {
    // Initialize the ears
    // There must always be  least two
    var ears = [Int]()

    // Get the status of the first two points
    // A point is an ear (and could be clipped) if the
    // line from point i-1 to i+1 is a diagonal of the polygon
    for i in 0...1 {
      if isDiagonal(first:i - 1, second:i + 1) {
        ears.append(i)
      }
    }

    // Quadrilaterals have a useful symmetry which means no more
    // diagonals need to be computed
    if polygonCount > 4 {
      // Get the remaining points
      for i in 2..<polygonCount {
        if isDiagonal(first: i - 1, second: i + 1) {
          ears.append(i)
        }
      }
    } else {
      // A quadrilateral is a special case due to symmetry
      // Just copy the first two entries (saves two calls to isDiagonal)
      for i in 2...3 {
        if ears.contains(i - 2) {
          ears.append(i)
        }
      }
    }

    // Scrambled polygon
    if ears.count < 2 {
      throw zoneError.initError("Zone \(Zone.name!) can't find two ears in the polygon \(zoneIndices)")
    }

    // Return the ear list
    return ears
  }
  
  // Returns the orientation of a triangle of polygon vertices i, j, k
  // Convex vertices will be positive
  func orientation(_ i: Int, _ j:Int, _ k:Int) -> Double {
    // Returns +1 for positive (left hand turn)
    //          0 for colinear
    //         -1 for negative (right hand turn)
    let a = 2 * zoneIndices[(i + polygonCount) % polygonCount]
    let b = 2 * zoneIndices[(j + polygonCount) % polygonCount]
    let c = 2 * zoneIndices[(k + polygonCount) % polygonCount]
    
    return orient(Triangulation.coords[a], Triangulation.coords[a + 1],
                        Triangulation.coords[b], Triangulation.coords[b + 1],
                        Triangulation.coords[c], Triangulation.coords[c + 1])
  }
  
  //  isDiagonal: Do two points lie on a diagonal of the polygon?
  func isDiagonal(first a: Int, second b: Int) -> Bool {
    // Two points a and b lie on a diagonal if the
    // line joining them lies within the polygon and the
    // points are mutually visible
    // return inCone(a, b) && inCone(b, a) && isVisible(a, b)
    return inCone(a, b) && isVisible(a, b)
  }
  
  //  inCone:  Is the diagonal line between points a & b within the cone
  //             formed by the sides of the polygon at the point ?
  func inCone(_ a:Int, _ b:Int) -> Bool {
    // Get the three points a - 1, a and a + 1
    let a0 = a - 1
    let a1 = a + 1
    
    // Test the orientation of the three triangles [a, a + 1, b], [a, b, a - 1] and [a, a - 1, a + 1]
    // Count how many are false
    var score = orientation(a, a1, b) > 0 ? 0 : 1
    
    // Test triangle [a, b, a0]
    if orientation(a, b, a0) <= 0 {
      score += 1
    }
    
    // Only need to test the last triangle [a, a0, a1] when count is unity
    // since if zero result must be < 2; when two result must be >= 2
    if score == 1 && orientation(a, a0, a1) <= 0 {
      return false
    }
    
    // The diagonal is within the cone if the score is less than two
    return score < 2
  }
  
  //  isVisible:  Is point (i) visible to point (j) (or vice-versa)
  //
  // For weak simple polygons this will fail because
  // there can be duplicate points
  func isVisible(_ i:Int, _ j:Int) -> Bool {
    // The vertices at these indices
    let u = zoneIndices[(i + polygonCount) % polygonCount], v = zoneIndices[(j + polygonCount) % polygonCount]
    
    // Get each edge
    for index in 0..<polygonCount {
      // An edge
      let f = [zoneIndices[index], zoneIndices[(index + 1) % polygonCount]]

      // Is the vertex at i or the vertex at j included in this edge?
      if !f.contains(u) && !f.contains(v) {
        // If the polygon side [index, index + 1] intersects the side [i, j] the vertices
        // are not inter-visible
        if edgeIntersectsEdge(u, v, f[0], f[1]) {
          return false
        }
      }
    }
    
    // The vertices are intervisible
    return true
  }

  // Do two edges s = [a, b] and t = [c, d] intersect?
  //
  func edgeIntersectsEdge(_ a:Int, _ b:Int, _ c:Int, _ d:Int) -> Bool {
    // These form four triangles
    //                    d
    //                    |
    //            a ------|-------- b
    //                    |
    //                    c
    //
    
    // See if the edges potentially intersect
    // Compute the orientations
    let ax = Triangulation.coords[2 * a], ay = Triangulation.coords[2 * a + 1]
    let bx = Triangulation.coords[2 * b], by = Triangulation.coords[2 * b + 1]
    let cx = Triangulation.coords[2 * c], cy = Triangulation.coords[2 * c + 1]
    let dx = Triangulation.coords[2 * d], dy = Triangulation.coords[2 * d + 1]
    
    let abcOrient = orient(ax, ay, bx, by, cx, cy)
    let abdOrient = orient(ax, ay, bx, by, dx, dy)
    let cdaOrient = orient(cx, cy, dx, dy, ax, ay)
    let cdbOrient = orient(cx, cy, dx, dy, bx, by)
    
    
    // If exactly one each of the pairs of triangles abc, abd & cda, cdb
    // are positive then the edges must intersect
    if ((abcOrient > 0) != (abdOrient > 0) && (cdaOrient > 0) != (cdbOrient > 0)) {
      return true
    }
    
    // The edges intersect if any combination of three points are colinear
    // and the third point lies 'between' the other two points
    return (abcOrient == 0 && isBetween(first: a, second: b, middle: c)) ||
      (abdOrient == 0 && isBetween(first: a, second: b, middle: d)) ||
      (cdaOrient == 0 && isBetween(first: c, second: d, middle: a)) ||
      (cdbOrient == 0 && isBetween(first: c, second: d, middle: b))
  }
} // End of zones

// Is the point p inside a closed loop of points?
func isInside(boundary indices:Array<Int>, size polygonCount:Int, query_x queryX:Double, query_y queryY:Double) -> Bool {
  // Count all the times a horizontal ray
  // heading to infinity from the query point
  // crosses the hull
  if polygonCount < 3 {
    return false
  }
  
  // Get each edge of the perimeter
  var crossings = 0
  for i in 0..<polygonCount {
    // Each edge is bounded by a pair of points p, q
    let j = (i + 1) % polygonCount
    let py = Triangulation.coords[2 * indices[i] + 1] - queryY,
        qy = Triangulation.coords[2 * indices[j] + 1] - queryY
    
    // When equal to zero need to include
    // one boundary point (p) and exclude
    // the other (q) otherwise
    // crossings will be double-counted
    if (py >= 0) != (qy > 0) {
      let px = Triangulation.coords[2 * indices[i]] - queryX,
          qx = Triangulation.coords[2 * indices[j]] - queryX
      
      // Sum area signs
      let area = orient(px, py, qx, qy, 0, 0)
      if area != 0 {
        crossings += area > 0 ? 1 : -1
      }
    } // End of if both same sign
  }
  
  // True if positive
  return crossings > 0
}

// Is a loop clockwise or anticlockwise?
func isClockwise(_ loop:Array<Int>) -> Bool {
  // Get the bottom right corner of the loop
  let a = loop.min { u, v in
    let ux = Triangulation.coords[2 * u], vx = Triangulation.coords[2 * v]

    if ux != vx { return ux < vx}
    let uy = Triangulation.coords[2 * u + 1], vy = Triangulation.coords[2 * v + 1]
    return uy < vy
  }
  
  // Which index is this
  let i = loop.firstIndex(of: a!)
  
  // Get the area of the triangle here
  let polygonCount = loop.count
  let b = loop[(i! + 1) % polygonCount]
  let c = loop[(i! + polygonCount - 1) % polygonCount]

  return orient(Triangulation.coords[2 * a!], Triangulation.coords[2 * a! + 1],
                Triangulation.coords[2 * b], Triangulation.coords[2 * b + 1],
                Triangulation.coords[2 * c], Triangulation.coords[2 * c + 1]) < 0
}

// Process the input data
// Can be either a zone file or a precomputed triangulation
// 
//
func triangulateZones(using storedData:StoredZones) throws {

  // Set the static variables
  if let c = storedData.conformingTo {
    Zone.hullConforming = (c  == "Voronoi")
  } else {
    Zone.hullConforming = true
  }
    
  //  Some control settings 
  if let c = storedData.control {
    Zone.control = Control(from: c)
  } else {
    Zone.control = Control()
  }
  
  // Need to establish a bounding box
  // We can do this here because
  // other points added at the instance level always
  // conform the the specified perimeter
  
  // Get the arcs (must be present)
  var arcs = storedData.arcs

  // And the properties (providing a default if none defined)
  let properties:Array<Properties> = storedData.properties ?? [Properties()] 
  
  // Get the transform
  if let transform = storedData.transform {
    // The arcs are quantized
    for i in 0..<arcs.count {
      // Current point (translated)
      var x:Double = 0, y:Double = 0
      
      // Get each subsequent point in the arc - scaling is done at a later stage
      for j in 0..<arcs[i].count {
        x += arcs[i][j][0]
        y += arcs[i][j][1]
        arcs[i][j][0] = x
        arcs[i][j][1] = y
      }
    }
    
    // The transform is the global scale
    Zone.scale = transform.scale!
    Zone.origin = transform.translate!
  }
 
  // The objects
  let objects:Dictionary<String, Array<StoredZones.GeometryCollection>> = storedData.objects
  var convexZoneList = Array<Zone>()
  for (key, geometryCollection) in objects {
    // Empty name? 
    if nil == Zone.name {
      if let n = storedData.name {
        Zone.name = n
      } else { 
        Zone.name = key 
      }
    }
 

    // Read each instance
    for g in geometryCollection {
      for z in g.geometries {
        convexZoneList.append(contentsOf:try! Zone(usingData: z, with: arcs).convexZones)
      }
    }
  }
  
  // Now triangulate
  
  // Make space for all the points
  Triangulation.pointCount = Triangulation.coords.count / 2
  Triangulation.triangulation = Triangulation(size: Triangulation.pointCount)
  
  if showMe != 0 {
    print("Pass 1 - Triangulate initial \(convexZoneList.count) zones")
  }
  
  for z in convexZoneList {
    // Add the zone to the triangulation
    do {
      let convexHullNext = try Triangulation.triangulation.addVertices(indices: z.zoneIndices, boundary: z.code)
      
      // Join convexZones together
      try Triangulation.triangulation.join(loop: convexHullNext)
    } catch {
      fatalError("Couldn't triangulate \(z)\n\(error)")
    }
  } // End of all convexzones
  
  // Need to compute the hull
  _ = Triangulation.triangulation.hullLoops()
  
  // Add circumvertices
  if Zone.hullConforming ?? true {
    try! Triangulation.triangulation.voronoiConforming()
  }
}


