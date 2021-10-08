//
//  ViewController.swift
//  Contact Book
//
//  Created by dato che on 12/24/20.
//

import UIKit
import CoreData

class ViewController: UIViewController {
    
    @IBOutlet var collectionView: UICollectionView!
    
    var dbContext = DBManager.shared.persistentContainer.viewContext
    
    lazy var people = [Person]()

    lazy var flowLayout: UICollectionViewFlowLayout = {
        let flowLayout = UICollectionViewFlowLayout()
        flowLayout.scrollDirection = .vertical
        return flowLayout
    }()
    
    
    override func viewDidLoad() {
        super.viewDidLoad()
        // Do any additional setup after loading the view.
        collectionView.delegate = self
        collectionView.dataSource = self
        
        collectionView.collectionViewLayout = flowLayout
        definesPresentationContext = true
        collectionView.register(
            UINib(nibName: "CollectionViewCell", bundle: nil),
            forCellWithReuseIdentifier: "PersonCell"
        )
        
        collectionView.addGestureRecognizer(
            UILongPressGestureRecognizer(
                target: self,
                action: #selector(handleLongPress(gesture:))
            )
        )
        
        fetchPeople()
    }
    
    func fetchPeople() {
          let request = Person.fetchRequest() as NSFetchRequest<Person>

        do {
            people = try dbContext.fetch(request)
            collectionView.reloadData()
        } catch {}
    }

    
    @objc func handleLongPress(gesture: UILongPressGestureRecognizer) {
        let location = gesture.location(in: collectionView)
        if let indexPath = collectionView.indexPathForItem(at: location) {
            deleteContact(indexPath)
        }
    }
    
    func deleteContact(_ indexPath: IndexPath){
        let person = people[indexPath.row]
        let name = person.name ?? ""
        let messageText = "Contact " + name + " will be deleted"
        let alert = UIAlertController(title: "Delete Contact", message: messageText, preferredStyle: .alert)
        
        alert.addAction(
            UIAlertAction(
                title: "Cancel",
                style: .cancel,
                handler: nil
            )
        )
        
        alert.addAction(
            UIAlertAction(
                title: "Delete",
                style: .destructive,
                handler: { [unowned self] _ in
                    dbContext.delete(person)
                    do {
                        try dbContext.save()
                        fetchPeople()
                    } catch {}
                }
            )
        )

        //ამ ადგილას წერს, რომ უკვე არის დაპრეზენტებულიო და თავიდან აპრეზენტებო და ვერ ვხვდები რა უნდა.
        //add-ის დროს არ აქვს მსგავსი პრობლემა
        present(alert, animated: true, completion: nil)
    }


    
    @IBAction func add() {
        let alert = UIAlertController(title: "Add Contact", message: nil, preferredStyle: .alert)
        
        var nameField: UITextField?
        var numberField: UITextField?
        
        alert.addTextField { textField in
            textField.placeholder = "Contact Name"
            nameField = textField
        }
        
        alert.addTextField { textField in
            textField.placeholder = "Contact Number"
            numberField = textField
        }
        
        alert.addAction(
            UIAlertAction(
                title: "Save",
                style: .default,
                handler: { [unowned self] _ in
                    guard let namefield = nameField, let numfield = numberField, let name = namefield.text, let number = numfield.text, !number.isEmpty, !name.isEmpty else { return }
                    
                    let person = Person(context: dbContext)
                    person.name = name
                    person.phone_number = number
                    
                    do {
                        try dbContext.save()
                        fetchPeople()
                    } catch {}
                }
            )
        )
        alert.addAction(
            UIAlertAction(
                title: "Cancel",
                style: .cancel,
                handler: nil
            )
        )
        present(alert, animated: true, completion: nil)
    }
}


extension ViewController: UICollectionViewDelegate {
    
}

extension ViewController: UICollectionViewDataSource {
    func numberOfSections(in collectionView: UICollectionView) -> Int {
        return 1
    }
    
    func collectionView(_ collectionView: UICollectionView, numberOfItemsInSection section: Int) -> Int{
        return people.count
    }
    
    func collectionView(_ collectionView: UICollectionView, cellForItemAt indexPath: IndexPath) -> UICollectionViewCell {
        let cell = collectionView.dequeueReusableCell(withReuseIdentifier: "PersonCell",
                                                      for: indexPath)
        if let personCell = cell as? CollectionViewCell{
            let person = people[indexPath.row]
            let wholeName = person.name
            let components = wholeName?.components(separatedBy: " ")
            let firstName = String(components?[0] ?? "")
        
            personCell.nameLabel.text = firstName
            let nameInitial = firstName[firstName.startIndex].uppercased()
            var surnameInitial = ""
            if components?.count ?? 0 > 1 {
                let surname = String(components?[1] ?? "")
                surnameInitial = surname[surname.startIndex].uppercased()
            }

        
            personCell.initialsLabel.text = nameInitial + surnameInitial

            personCell.numberLabel.text = person.phone_number

        }
        return cell
    }

    
}

extension ViewController: UICollectionViewDelegateFlowLayout {
    
    func collectionView(
        _ collectionView: UICollectionView,
        layout collectionViewLayout: UICollectionViewLayout,
        insetForSectionAt section: Int
    ) -> UIEdgeInsets {
        return UIEdgeInsets(top: Constants.spacing, left: Constants.spacing, bottom: Constants.spacing, right: Constants.spacing)
    }
    
    func collectionView(
        _ collectionView: UICollectionView,
        layout collectionViewLayout: UICollectionViewLayout,
        sizeForItemAt indexPath: IndexPath
    ) -> CGSize {
        let spareWidth = collectionView.frame.width - (2 * Constants.spacing) - ((Constants.itemCountInLine - 1) * Constants.spacing)
        return CGSize(width: spareWidth / Constants.itemCountInLine, height: collectionView.frame.height / 5)
    }
    
    func collectionView(
        _ collectionView: UICollectionView,
        layout collectionViewLayout: UICollectionViewLayout,
        minimumInteritemSpacingForSectionAt section: Int
    ) -> CGFloat {
        return Constants.spacing
    }
    
    func collectionView(
        _ collectionView: UICollectionView,
        layout collectionViewLayout: UICollectionViewLayout,
        minimumLineSpacingForSectionAt section: Int
    ) -> CGFloat {
        return Constants.heightSpacing
    }
}

extension ViewController {
    
    struct Constants {
        static let itemCountInLine: CGFloat = 3
        static let spacing: CGFloat = UIScreen.main.bounds.width * 0.04
        static let heightSpacing: CGFloat = UIScreen.main.bounds.height * 0.05
    }
}
