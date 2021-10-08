//
//  Person+CoreDataProperties.swift
//  Contact Book
//
//  Created by dato che on 12/24/20.
//
//

import Foundation
import CoreData


extension Person {

    @nonobjc public class func fetchRequest() -> NSFetchRequest<Person> {
        return NSFetchRequest<Person>(entityName: "Person")
    }

    @NSManaged public var phone_number: String?
    @NSManaged public var name: String?

}

extension Person : Identifiable {

}
