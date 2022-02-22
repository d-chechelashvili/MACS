package ge.dchechelashvili.notesapp.details

import android.view.KeyEvent
import android.view.LayoutInflater
import android.view.View
import android.view.ViewGroup
import android.widget.CheckBox
import android.widget.EditText
import android.widget.ImageView
import androidx.recyclerview.widget.RecyclerView
import ge.dchechelashvili.notesapp.R
import ge.dchechelashvili.notesapp.data.entity.Item


class CheckListAdapter(val listListener: CheckListListener) :
    RecyclerView.Adapter<CheckListAdapter.CheckListViewHolder>() {

    var items = mutableListOf<Item>()

    override fun onCreateViewHolder(
        parent: ViewGroup,
        viewType: Int
    ): CheckListAdapter.CheckListViewHolder {
        return CheckListViewHolder(
            LayoutInflater.from(parent.context).inflate(R.layout.checklist_item, parent, false)
        )
    }

    override fun onBindViewHolder(holder: CheckListAdapter.CheckListViewHolder, position: Int) {
        val item = items[position]
        holder.bindItem(item)
        holder.checkBox.setOnClickListener {
            listListener.onCheckBoxClicked(item, position)
        }
        holder.contentText.setOnKeyListener(View.OnKeyListener { v, keyCode, event ->
            if (keyCode == KeyEvent.KEYCODE_ENTER && event.action == KeyEvent.ACTION_UP) {
                val newItem = Item(holder.contentText.text.toString(), item.noteId, item.checked)
                items.set(position, newItem)
                return@OnKeyListener true
            }
            false
        })
        holder.removeButton.setOnClickListener{
            listListener.onRemoveButtonClicked(item, position)
        }
    }

    override fun getItemCount(): Int {
        return items.size
    }

    inner class CheckListViewHolder(view: View) : RecyclerView.ViewHolder(view) {
        fun bindItem(item: Item) {
            checkBox.isChecked = item.checked
            contentText.setText(item.content)
            if (checkBox.isChecked){
                contentText.focusable = View.NOT_FOCUSABLE
                removeButton.visibility = View.GONE
            }else{
                contentText.focusable = View.FOCUSABLE
                removeButton.visibility = View.VISIBLE
            }
        }

        val checkBox: CheckBox = view.findViewById(R.id.itemCheckBox)
        val contentText: EditText = view.findViewById(R.id.editTextContent)
        val removeButton: ImageView = view.findViewById(R.id.removeButton)

    }


}