//  VoroBox
//
//  StoredZones.swift
//
//  Created by Rob Whitehurst on 13/6/20.
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

import Foundation

// The Basic structure for JSON IO
struct StoredZones : Codable {
  var zones: [Zone]
  var points: Dictionary<String, Array<Double>>?
  var coordinates: Array<Array<Double>>?

  var name: String
  var conformingTo:String?
  
  // global properties
  var properties: Properties?
  var id:Int?
  var iteration:Int?

  var origin: Array<Double>?
  var scale: Double?
  var order: String? = "anticlockwise"
  
  struct Zone: Codable {
    var name: String?

    var origin: Array<Double>?
    var scale: Double?
    var order: String? = "anticlockwise"
    
    // Standard description
    var boundary:Array<String>?
    var holes:Array<Array<String>>?

    // Alternative geojson-like description
    var polygon:Array<Array<Array<Double>>>?
    var multipolygon:Array<Array<Array<Array<Double>>>>?

    // Properties associated with the zone
    var properties: Properties?
  }
  

  
  struct Properties: Codable {
    // Zone properties
    var randomDensity:Double?
    var id:Int?
    var distinct:Bool?

    // debugging control
    var showMe:Int?
  }
}


