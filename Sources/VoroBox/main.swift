//
//  main.swift
//  VoroBox
//
//  Created by Rob Whitehurst on 2/7/20.
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
import Files
import Yams

let ZoneFile   = "zoneData.yml"

// Read input data using JSON
do {
  try triangulateZones(using: load(ZoneFile))
} catch {
  print ("Zone creation failed with error: \(error)")
}
// Write output data
var filename = "triangleData.json"
do {
  let iteration:Int = Zone.iteration + 1
  let folder = try Folder(path: OutputFolder)
  let conforming = Zone.hullConforming ?? true
  
  // File names - First the filled cells
  var filename = "cells_" + Zone.name
  if iteration > 1 {
    filename += "\(iteration)"
  }
  var file = try folder.createFile(named: filename + ".vtk")
  if conforming  {
    let x = Triangulation.triangulation.toVoronoi()
    try file.write(Triangulation.triangulation.vtkVoronoi(x))
    
    file = try folder.createFile(named: Zone.name + ".yml")
    try file.write(Triangulation.triangulation.toZone())
  } else {
    // Delaunay
    try file.write(Triangulation.triangulation.vtkDelaunay())
  }
  
  // And a vtk file showing the edges
  filename = "edges_" + Zone.name
  if iteration > 1 {
    filename += "\(iteration)"
  }
  file = try! folder.createFile(named: filename + ".vtk")
  try file.write(Triangulation.triangulation.vtkEdges(label: Zone.name, type: "Voronoi"))
  if !conforming {
    file = try! folder.createFile(named: filename + "_triangulation.vtk")
    try file.write(Triangulation.triangulation.vtkEdges(label: Zone.name, type: "Delaunay"))
  }
} catch {
  fatalError("Couldn't write \(filename)\n\(error)")
}


//let triangulationData: [Triangulation] = zoneData.map({$0.triangulation ?? Triangulation()})
func load<T: Decodable>(_ filename: String) -> T {
  do {
    // Locate the zone file
    let file = try Folder(path: ZoneFolder).file(named: filename)

    // read the data
    let data = try file.read()

    // Unpack the json  and
    // return the decoded data
    
    if ZoneFile == "zoneData.yml" {
      return try YAMLDecoder().decode(T.self, from: data)
    }
    return try JSONDecoder().decode(T.self, from: data)
  } catch {
    fatalError("Couldn't parse \(filename) as \(T.self):\n\(error)")
  }
}
