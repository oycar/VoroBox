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
  var randomDensity:Double?
  var id:Int?
  var distinct:Bool?
  var showMe:Int?
}


extension Properties {
  init(with density:Double? = nil, _ i:Int? = nil, _ flag:Bool? = true, _ show:Int? = nil) {
    randomDensity = density
    id = i
    distinct = flag
    showMe = show
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
  static var name:String = "Zone Name"
  
  // This should be a Trisngulation Boolean
  static var hullConforming: Bool?

  // A counter - start at FirstZone
  static var zoneID = FirstZone
  
  // Iteration
  static var iteration:Int = 0
  
  // Labelled Points
  static var labelledPoints = Dictionary<String, Int>()
  
  // The optional instance name
  var name: String?

  // The polygon bounding the zone
  // along with its size
  var polygon = Array<String>()
  var polygonCount = 0

  // Convex zones
  var convexZones = Array<Zone>()
  
  // Indices are associated with each instance
  var zoneIndices = Array<Int>()
    
  // A code
  var code: Int = 1

  // Control of point insertion
  static var globalProperties:Properties?
  var properties: Properties?
}

// The actual zone structure ...
extension Zone {
  // Read in a zone defined in a data file
  init(usingData stored: StoredZones.Zone) throws {
    // Name this zone
    name = stored.name
    
    // when showMesh is set lots of step by step snapshots of the mesh are drawn

    
    // Global properties can give defaults
    var distinct = Zone.hullConforming ?? true
    if let g = Zone.globalProperties {
      properties = Properties(with: g.randomDensity, g.id, g.distinct, g.showMe)
      if nil != g.distinct { distinct = g.distinct! }
      if nil != g.showMe { showMe = g.showMe! }
      
      if let p = stored.properties {
        properties = Properties(with: p.randomDensity ?? g.randomDensity, p.id ?? g.id)
        if nil != p.distinct { distinct = p.distinct! }
      }
    } else if let p = stored.properties {
      // No default
      properties = Properties(with: p.randomDensity, p.id)
      if nil != p.distinct { distinct = p.distinct! }
      if nil != p.showMe { showMe = p.showMe! }
    }
    
    
    // Set the zone code
    code = 2 * Zone.zoneID

    // Update the code depending on if this zone is distinct or not
    if distinct { code += 1 } // Code is odd for conforming edges
    
    // Set code - when even the zone is Voronoi conforming (the default)
    Zone.zoneID += 1
    
    // Initially we are reading the stored data, which can define
    // a reflex polygon - so we need to create a collection
    // of convex zones
    
    // We are given the co-ordinates
    // Either as labelled points in a list of strings
    // Need to find the points that correspond to these
    if var list = stored.list {
      
      // Check ordering
      if stored.order == "clockwise" { list.reverse() }
      
      // Add the indices
      for label in list {
        let i = Zone.labelledPoints[label]!
        zoneIndices.append(i)
      }
      
      // The polygon order
      polygonCount = zoneIndices.count
    } else if var polygon = stored.polygon {
      // Or as explicit polygon of [x,y] points
      // Check ordering
      if stored.order == "clockwise" { polygon.reverse() }

      // In this case
      // Zones can be referred to a specific origin
      var origin  = Array<Double>([0, 0])
      if let o = stored.origin {
        origin = o
      }
      
      // And a scale factor can be applied
      var scale = 1.0
      if let s = stored.scale {
        scale = s
      }
      
      // The points are mapped to vertices
      roundPolygon: for k in 0..<polygon.count {
        var p = polygon[k]
        for j in 0...1 {
          p[j] = origin[j] + scale * p[j]
        }
        
        for i in 0..<Triangulation.pointCount {
          let x = dist(Triangulation.coords[2 * i], Triangulation.coords[2 * i + 1], p[0], p[1])
          
          // Is this a matching point?
          if x < squaredThreshold {
            zoneIndices.append(i)
            continue roundPolygon
          }
        }
        
        // A new point if none found - input points are treated as members of a triangulation
        Triangulation.coords.append(contentsOf: p)
        
        // Use the zone properties
        Triangulation.code.append(NoZoneCode)
        
        zoneIndices.append(Triangulation.pointCount)
        Triangulation.pointCount += 1
      }
      polygonCount = zoneIndices.count

    }
    
    if !zoneIndices.isEmpty {
    
      // Now need to make a list of convex zones
      do {
        try clipZone()

      } catch {
        fatalError("Couldn't triangulate \(self)\n\(error)")
      }
    } // Clipped the zone into convex polygons
    
    // The final convex zone is simply self
    convexZones.append(self)
  } // init
 
  // Some specialized inits
  init(clip zone: Zone, ear i:Int, order n:Int, label optionalName:String? = nil) {
    // This handles the case of convex zones
    // This generates a triangular zone
    name = optionalName
    polygonCount = n // Its an n-gon - can be zero
    
    // clipped zones inherit parent code
    code = zone.code
    
    // Needed for initialization
    properties = zone.properties
    
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
  
  mutating func addPointsToZone(_ pointList: Array<Int>) -> Array<Int> {
    // Now we can allocate the points (if any) that are associated with each convex zone
    var newList = Array<Int>()
    
    // Add defined points
    for i in pointList {
      // Do we have an unattached specified point?
      let x = Triangulation.coords[2 * i], y = Triangulation.coords[2 * i + 1]
      
      // Is this inside this zone?
      if polygonCount == 0 || isInside(query_x: x, query_y: y) {
        // Yes this point is inside this zone
        zoneIndices.append(i)
        
        // Can only be in one zone
      } else {
        newList.append(i)
      }
    }
    
    // All done if no perimeter or bounding box
    if polygonCount < 2 { return newList}
    
    // Add points defined by a perimter
    var i = 1
    var minX = Triangulation.coords[2 * zoneIndices[0]]
    var minY = Triangulation.coords[2 * zoneIndices[0] + 1]
    var maxX = minX
    var maxY = minY
    
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
    
    // Now add random points
    // Update zone count
    Triangulation.pointCount = Triangulation.coords.count / 2
    if let density:Double = properties!.randomDensity {
      // Number of points is density times box area
      var r:Double = density * (maxX - minX) * (maxY - minY)
      r.round(.up)
      
      // Add the points
      for _ in 0..<Int(r) {
        // Add points inside  the zone
        let x = Double.random(in: minX...maxX), y = Double.random(in: minY...maxY)
        
        // Append the point
        if isInside(query_x: x, query_y: y) {
          // New point - set coordinates
          Triangulation.coords.append(x)
          Triangulation.coords.append(y)
          
          // This needs the zone code
          Triangulation.code.append(code)
          
          // Set zone index
          zoneIndices.append(Triangulation.pointCount)
          Triangulation.pointCount += 1
        }
      }
    } // Add random points
    
    return newList
  }
  
  // Is the point p inside the current zone?
  func isInside(query_x queryX:Double, query_y queryY:Double) -> Bool {
    // Check each edge in the polygon hull
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
      let py = Triangulation.coords[2 * zoneIndices[i] + 1] - queryY,
      qy = Triangulation.coords[2 * zoneIndices[j] + 1] - queryY
      
      // When equal to zero need to include
      // one boundary point (p) and exclude
      // the other (q) otherwise
      // crossings will be double-counted
      if (py >= 0) != (qy > 0) {
        let px = Triangulation.coords[2 * zoneIndices[i]] - queryX,
        qx = Triangulation.coords[2 * zoneIndices[j]] - queryX
        
        // Sum area signs
        let area = orientIfSure(px, py, qx, qy, 0, 0)
        if area != 0 {
          crossings += area > 0 ? 1 : -1
        }
      } // End of if both same sign
    }
    
    // True if positive
    return crossings > 0
  }
  
  // Extra geometric primitives for zones
  
  
  // Clip the zone
  mutating func clipZone() throws {
    // Anything less than quadrilateral is convex
    if polygonCount < 4 { return }
    
    // get the initial ear list
    var ears:Array<Int> = try listEars()
    
    // Determine which ears should be clipped
    // Finished when all are ears
    // Keep track of the last clipped zone
    var convexZone:Zone? = nil
    untilConvex: while polygonCount > 3 && ears.count < polygonCount {
      // For testing turn off reverible clip
      var reversibleClip = true
      
      // Get the ear at the end of the list
      let i = ears.last!
      
      // Clip corner i
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
    throw zoneError.initError("Zone \(Zone.name) can't find matching vertices to join two zones")
  }
  
  //
  // ears :  Compute the ear status of all points
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
      throw zoneError.initError("Zone \(Zone.name) can't find two ears in the polygon \(zoneIndices)")
    }
    
    // Return the ear list
    return ears
  }
  
  // Returns the orientation of a triangle of polygon vertices i, j, k
  // Convex vertices will be positive
//  func orientation(_ i: Int, _ j:Int, _ k:Int) -> Double {
//    // Returns +1 for positive (left hand turn)
//    //          0 for colinear
//    //         -1 for negative (right hand turn)
//    let a = 2 * zoneIndices[(i + polygonCount) % polygonCount]
//    let b = 2 * zoneIndices[(j + polygonCount) % polygonCount]
//    let c = 2 * zoneIndices[(k + polygonCount) % polygonCount]
//    return orientIfSure(Triangulation.coords[a], Triangulation.coords[a + 1],
//                        Triangulation.coords[b], Triangulation.coords[b + 1],
//                        Triangulation.coords[c], Triangulation.coords[c + 1])
//  }
  func orientation(_ i: Int, _ j:Int, _ k:Int) -> Double {
    // Returns +1 for positive (left hand turn)
    //          0 for colinear
    //         -1 for negative (right hand turn)
    let a = 2 * zoneIndices[(i + polygonCount) % polygonCount]
    let b = 2 * zoneIndices[(j + polygonCount) % polygonCount]
    let c = 2 * zoneIndices[(k + polygonCount) % polygonCount]
    
    return orientIfSure(Triangulation.coords[a], Triangulation.coords[a + 1],
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
  func isVisible(_ i:Int, _ j:Int) -> Bool {
    let ij = [(i + polygonCount) % polygonCount, (j + polygonCount) % polygonCount]
    for index in 0..<polygonCount {
      // An edge
      let e = [index, (index + 1) % polygonCount]
      
      // Is i or j included in this edge?
      if !e.contains(ij[0]) && !e.contains(ij[1]) {
        // If the polygon side [index, index + 1] intersects the side [i, j] the vertices
        // are not inter-visible
        if edgeIntersectsEdge(firstEdge: ij, secondEdge: e) {
          return false
        }
      }
    }
    
    // The vertices are intervisible
    return true
  }
  
  // Do two edges s = [a, b] and t = [c, d] intersect?
  //
  func edgeIntersectsEdge(firstEdge s: Array<Int>, secondEdge t: Array<Int>) -> Bool {
    // Obtain the four points a, b, c and d located on each side
    let a = s[0]
    let b = s[1]
    let c = t[0]
    let d = t[1]
    
    // These form four triangles
    //                    d
    //                    |
    //            a ------|-------- b
    //                    |
    //                    c
    //
    
    
    // See if the edges potentially intersect
    // Compute the orientations
    let abcOrient = orientation(a, b, c)
    let abdOrient = orientation(a, b, d)
    let cdaOrient = orientation(c, d, a)
    let cdbOrient = orientation(c, d, b)
    
    // Assume the sides do not intersect
    var intersect = false
    
    // If exactly one each of the pairs of triangles abc, abd & cda, cdb
    // are positive then the edges must intersect
    if ((abcOrient > 0) != (abdOrient > 0) && (cdaOrient > 0) != (cdbOrient > 0)) {
      intersect = true
    }
    
    // The edges intersect if any combination of three points are colinear
    // and the third point lies 'between' the other two points
    // This never gets called?
    if !intersect { // If intersection is known nothing to do
      // More efficient to precompute modified labels
      let p = zoneIndices[(a + polygonCount) % polygonCount]
      let q = zoneIndices[(b + polygonCount) % polygonCount]
      let r = zoneIndices[(c + polygonCount) % polygonCount]
      let s = zoneIndices[(d + polygonCount) % polygonCount]
      
      // Just compute the intersectioq
      intersect = (abcOrient == 0 && isBetween(first: p, second: q, middle: r)) ||
        (abdOrient == 0 && isBetween(first: p, second: q, middle: s)) ||
        (cdaOrient == 0 && isBetween(first: r, second: s, middle: p)) ||
        (cdbOrient == 0 && isBetween(first: r, second: s, middle: q))
    }
    
    // Return the intersection status
    return intersect
  }

  
    
//
//  // For screening input files
//  func checkPolygon() -> Bool {
//    // Get each edge index i
//    for i in 0..<zoneIndices.count-1 {
//      let jmax = i == 0 ? zoneIndices.count - 1 : zoneIndices.count
//      for j in i+2..<jmax {
//        let a = [i, (i + polygonCount) % polygonCount]
//        let b = [j, (j + polygonCount) % polygonCount]
//
//        // Do the two edges i & j intersect?
//        if edgeIntersectsEdge(firstEdge: a, secondEdge: b) {
//          // Yes
//          return false
//        }
//      }
//    }
//    return true
//  }
} // End of zones

// Process the input data
// Can be either a zone file or a precomputed triangulation
// 
//
func triangulateZones(using storedData:StoredZones) throws {

  // Set the static variables
  // First the name
  Zone.name = storedData.name
  if let c = storedData.conformingTo {
    Zone.hullConforming = (c  == "Voronoi")
  } else {
    Zone.hullConforming = true
  }
  
  // Update zone ID
  if let storedID = storedData.id {
    if storedID > Zone.zoneID  { Zone.zoneID = storedID }
  }
  
  // Update iteration
  Zone.iteration = storedData.iteration ?? 0
  
  // Global properties
  if let p = storedData.properties {
    Zone.globalProperties = Properties(with: p.randomDensity, p.id, p.distinct, p.showMe)
  } else {
    Zone.globalProperties = Properties(with: nil, 0, true)
  }

  // Need to establish a bounding box
  // We can do this here because
  // other points added at the instance level always
  // conform the the specified perimeter
  var minX = boundingBox[0][0]
  var minY = boundingBox[0][1]
  var maxX = boundingBox[1][0]
  var maxY = boundingBox[1][1]
  
  // Now any labelled points (which make up perimeters etc)
  if let storedPoints = storedData.points {
    for (label, point) in storedPoints {
      let x = point[0], y = point[1]
      
      // Add to the list of points in the zone - always
      Triangulation.coords.append(x)
      Triangulation.coords.append(y)
      Triangulation.code.append(NoZoneCode)
      
      // Need to tag these points
      Zone.labelledPoints[label] = Triangulation.pointCount
      Triangulation.pointCount += 1

      // Compute bounds
      if (x < minX)  {minX = x}
      if (y < minY)  {minY = y}
      if (x > maxX)  {maxX = x}
      if (y > maxY)  {maxY = y}
    }
  }
  
  // Read each instance
  var convexZoneList = Array<Zone>()
  for (_, z) in storedData.zones.enumerated() {
    convexZoneList.append(contentsOf:try! Zone(usingData: z).convexZones)
  }
  
  // Note the current zone code
  let zCode = convexZoneList.last!.code
  
  // Scan the added points
  for i in 0..<Triangulation.pointCount {
    let x = Triangulation.coords[2 * i], y = Triangulation.coords[2 * i + 1]
        
    // Compute bounds
    if (x < minX)  {minX = x}
    if (y < minY)  {minY = y}
    if (x > maxX)  {maxX = x}
    if (y > maxY)  {maxY = y}
  }
  
  // Now the explicitly specified points
  var specifiedPoints = Array<Int>()
  if let storedPoints = storedData.coordinates {
    for (_, point) in storedPoints.enumerated() {
      let x = point[0], y = point[1]

      // Add it to the coordinate list
      Triangulation.coords.append(x)
      Triangulation.coords.append(y)

      // These are given a zone code
      Triangulation.code.append(zCode)
      
      specifiedPoints.append(Triangulation.pointCount)
      Triangulation.pointCount += 1

      // Compute bounds - make no assumption about labelled points
      if (x < minX)  {minX = x}
      if (y < minY)  {minY = y}
      if (x > maxX)  {maxX = x}
      if (y > maxY)  {maxY = y}
    }
  }
  
  // Update the overall bounding box
  if (boundingBox[0][0] > minX) { boundingBox[0][0] = minX }
  if (boundingBox[0][1] > minY) { boundingBox[0][1] = minY }
  if (boundingBox[1][0] < maxX) { boundingBox[1][0] = maxX }
  if (boundingBox[1][1] < maxY) { boundingBox[1][1] = maxY }
    
  // Set an initial triangulation - for the perimeters only at this point
  // These are all done in separate blocks deliberately
  // to minimize point shuffling in initializing the zones
  
  // Add the points to all the convex zones
  for (i, _) in convexZoneList.enumerated() {
    // Add the points
    specifiedPoints = convexZoneList[i].addPointsToZone(specifiedPoints)
  }

  // Make space for all the points
  Triangulation.pointCount = Triangulation.coords.count / 2
  Triangulation.triangulation = Triangulation(size: Triangulation.pointCount)

  // For reproducible runs output the initial zone file
  if showMe > 1 {
    // A reproducible zone file for debugging
    Triangulation.triangulation.showZone(zones: convexZoneList)
  }
  
  if showMe > 0 {
    print("Phase 0 - Triangulate initial \(convexZoneList.count) zones")
  }
  
  // Now triangulate
  for z in convexZoneList {
    // Add the zone to the triangulation
    do {
      let convexHullNext = try Triangulation.triangulation.addVertices(indices: z.zoneIndices, boundary: z.code)

      // Join convexZones together
      try Triangulation.triangulation.join(loop: convexHullNext)
      
      if showMe > 0 {
        print("Joined convex zone \(z)")
        Triangulation.triangulation.showVoronoi(label: "iz_", type: "Voronoi")
      }
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


